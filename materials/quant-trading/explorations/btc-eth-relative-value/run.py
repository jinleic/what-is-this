#!/usr/bin/env python3
"""BTC/ETH relative-value exploration engine (wave one).

Deterministic, stdlib-only. Windows:
  development 2020-08-01 .. 2023-12-31
  validation  2024-01-01 .. 2024-12-31
Fail-closed: --end beyond 2024-12-31 is refused unless QUANT_HOLDOUT_UNLOCK=1
is set in the environment (wave one never sets it).

Conventions (frozen before any result was inspected):
  * Signal inputs: hourly bars with close_time <= 00:00 UTC of day D, i.e.
    the last usable input is the daily close of D-1. Nothing after.
  * Fill: first trade print of day D (open of the 00:00 hourly bar).
  * Daily strategy return R_D: NAV just before the day D+1 trade (after the
    funding event that settles at 00:00 D+1 and is charged to day D's
    position, i.e. to the holder of the settled interval) divided by NAV
    just after the day D trade, minus 1. The final day of a window uses the
    truncated form NAV(close D)/NAV(fill D) - 1 so that no observation
    beyond --end is ever touched.
  * Funding: event with settlement timestamp T is charged to the position
    that held the settled interval. Events within 60 s of a UTC midnight
    settle an interval that ended at that midnight and are charged to the
    previous day's position; all other events to their own UTC date.
    Cashflow = -units_perp * rate * mark_price. Long pays positive funding.
  * Capital: fully collateralized. cash >= sum(|perp notional|) is asserted
    after every trade. No borrowing, no interest on cash, no leverage:
      long-winner/cash: 1.0 spot winner, rest cash.
      dollar-neutral:   0.5 spot long winner + 0.5 perp short loser
                        (gross 1.0, net 0.0, cash 1.0).
  * Costs: 12 bps one-way per unit notional on spot legs, 7 bps on
    perpetual legs (contract.json defaults); --stress multiplies both.
  * Benchmarks, (re)initialized at the start of each evaluated window:
      ew_spot_bah   equal-weight BTC/ETH spot buy-and-hold, 12 bps entry per
                    leg once, drift weights thereafter.
      static_dn     long 0.5 BTC spot + short 0.5 ETH perpetual, entered
                    once, never rebalanced; funding accrues on the short
                    leg, entry costs 12+7 bps, quantity fixed, NAV marks daily.
"""
import argparse
import csv
import io
import json
import math
import os
import sys
import zipfile
from datetime import date, datetime, timedelta, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "data", "raw", "binance"))
MS_MONTH = 30 * 24 * 3600 * 1000  # loose bucket only for month file iteration
WAVE_ONE_LAST_DAY = date(2024, 12, 31)
DEV_END = date(2023, 12, 31)

SPOT_BPS = 12.0 / 1e4
PERP_BPS = 7.0 / 1e4
ANN = 365.0  # crypto trades daily


# ---------------------------------------------------------------- data layer
def month_iter(start, end):
    y, m = start.year, start.month
    while (y, m) <= (end.year, end.month):
        yield y, m
        m += 1
        if m == 13:
            y, m = y + 1, 1


def read_zip_csv(path):
    with zipfile.ZipFile(path) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as fh:
            text = io.TextIOWrapper(fh, "utf-8")
            for row in csv.reader(text):
                if row:
                    yield row


def norm_ms(ts):
    ts = int(ts)
    if ts > 10 ** 14:  # microseconds variant (some spot months)
        ts //= 1000
    return ts


def load_hourly(symbol, kind, start, end, log):
    """Return {utc_ms_of_hour_open: (open,close)} hourly close map.

    kind: 'spot' or 'perp'. Only bars with open_time in [start, end] days
    are kept; nothing beyond --end is loaded or retained.
    """
    subdir = "spot-1h" if kind == "spot" else "um-perp-1h"
    out = {}
    n_rows = 0
    for y, m in month_iter(start, end):
        fname = "%s-1h-%04d-%02d.zip" % (symbol, y, m)
        path = os.path.join(DATA_ROOT, subdir, symbol, fname)
        if not os.path.exists(path):
            log.append("missing_kline %s %s" % (kind, fname))
            continue
        first = True
        for row in read_zip_csv(path):
            if first and not _is_int(row[0]):
                first = False  # header row (um-perp variant)
                continue
            first = False
            ts = norm_ms(row[0])
            d = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc)
            if d.date() >= start and d.date() <= end:
                out[ts] = (float(row[1]), float(row[4]))
                n_rows += 1
    log.append("loaded %s %s bars=%d" % (symbol, kind, n_rows))
    return out


def load_funding(symbol, start, end, log):
    """Return list of (settle_ms, rate) events with settle in [start, end+1d)."""
    out = []
    end_incl = end + timedelta(days=1)  # boundary event charged to last day
    for y, m in month_iter(start, end):
        fname = "%s-fundingRate-%04d-%02d.zip" % (symbol, y, m)
        path = os.path.join(DATA_ROOT, "funding-rate", symbol, fname)
        if not os.path.exists(path):
            log.append("missing_funding %s" % fname)
            continue
        for row in read_zip_csv(path):
            if not _is_int(row[0]):
                continue  # header
            ts = norm_ms(row[0])
            rate_col = row[2] if len(row) >= 3 else row[-1]
            try:
                rate = float(rate_col)
            except ValueError:
                continue
            dt = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc)
            if start <= dt.date() < end_incl:
                out.append((ts, rate))
    out.sort()
    log.append("loaded funding %s events=%d" % (symbol, len(out)))
    return out


def _is_int(tok):
    try:
        int(tok)
        return True
    except ValueError:
        return False


def daily_bars(hourly, log, label):
    """Aggregate an hourly map to {iso_date: (open, close)} using UTC days."""
    days = {}
    for ts in sorted(hourly):
        o, c = hourly[ts]
        d = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc).date().isoformat()
        if d not in days:
            days[d] = [o, c, ts]
        else:
            b = days[d]
            b[1] = c
            b[2] = ts
    drop = 0
    for d, b in list(days.items()):
        y, m, dd = (int(x) for x in d.split("-"))
        expected = (date(y, m, dd) - date(2020, 8, 1)).days * 24 + 24
        if b[2] == 0:
            drop += 1
    # completeness: report days with < 24 bars but keep them (close still valid)
    log.append("daily %s days=%d" % (label, len(days)))
    return {k: (v[0], v[1]) for k, v in sorted(days.items())}


def funding_by_day(events, log):
    """Map funding events to the UTC day whose position held the settled
    interval. Events within 60 s after a UTC midnight settle an interval that
    ended at that midnight: charged to the previous day's position."""
    out = {}
    for ts, rate in events:
        dt = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc)
        d = dt.date()
        if dt.hour == 0 and dt.minute == 0 and dt.second == 0 and ts % 60000 < 60000:
            d = d - timedelta(days=1)  # boundary settlement -> prior holder
        key = d.isoformat()
        out[key] = out.get(key, 0.0) + rate
    log.append("funding days mapped=%d" % len(out))
    return out


def perp_mark_lookup(permph, settle_ms):
    """Close of the hourly perp bar ending at settle time; fallback to the
    nearest earlier bar within 24h."""
    hour = (settle_ms // 3600000) * 3600000
    bar_open = hour - 3600000
    for k in range(24):
        got = permph.get(bar_open - k * 3600000)
        if got is not None:
            return got[1]
    return None


# ---------------------------------------------------------------- simulation
class Ledger(object):
    """Quantity-space ledger over a single symbol set.

    units: {'btc_spot': q, 'eth_spot': q, 'btc_perp': q, 'eth_perp': q}
    Fully collateralized: cash >= |perp notional| enforced after each trade.
    """

    def __init__(self, nav0=1.0):
        self.cash = nav0
        self.units = {"btc_spot": 0.0, "eth_spot": 0.0, "btc_perp": 0.0,
                      "eth_perp": 0.0}

    def mark(self, prices):
        """prices: dict slot->price (spot legs at spot prices, perp at perp)."""
        return self.cash + sum(self.units[s] * prices[s] for s in self.units)

    def perp_notional(self, prices):
        return (abs(self.units["btc_perp"] * prices["btc_perp"])
                + abs(self.units["eth_perp"] * prices["eth_perp"]))

    def trade_to(self, target_notional, prices, stress):
        """target_notional maps slot -> USD notional value. Quantities are
        derived as notional/price; shorting a perp produces negative units.
        Collateral rule: short perp notional must stay covered by cash (the
        fully collateralized constraint). Fees may push cash a hair negative
        on spot-only rebalances; that deficit is part of NAV (mark() includes
        it) and is not borrowing against open short risk."""
        spot_traded = 0.0
        perp_traded = 0.0
        for slot in sorted(self.units):
            p = prices[slot]
            if p <= 0:
                continue
            tgt_q = target_notional.get(slot, 0.0) / p
            delta_notional = (tgt_q - self.units[slot]) * p
            if abs(delta_notional) < 1e-12:
                continue
            self.units[slot] += delta_notional / p
            if slot.endswith("_spot"):
                self.cash -= delta_notional
                spot_traded += abs(delta_notional)
            else:
                if delta_notional < 0:
                    self.cash -= delta_notional  # short sale adds proceeds
                else:
                    self.cash -= delta_notional  # closing short pays out
                perp_traded += abs(delta_notional)
        cost = (spot_traded * SPOT_BPS + perp_traded * PERP_BPS) * stress
        self.cash -= cost
        turnover = spot_traded + perp_traded
        if self.perp_notional(prices) > 1e-12 and self.cash + 1e-9 < self.perp_notional(prices):
            raise RuntimeError("collateral breach: cash %.6f < short perp %.6f"
                               % (self.cash, self.perp_notional(prices)))
        if self.mark(prices) <= 0:
            raise RuntimeError("NAV wiped out")
        return turnover, cost

    def funding_event(self, slot, rate, mark_price):
        self.cash += -self.units[slot] * rate * mark_price


def weights_for_rule(mech, form, mode, r_btc, r_eth, ratio_z, prev_state):
    """Target unit weights for one daily rebalance. Per-leg caps 0.5 keep
    gross <= 1.0 (fully collateralized, no leverage).

    mech 'rm': relative momentum. Long the trailing-L winner.
      form 'cash': 1.0 winner spot, rest cash.
      form 'dn':   0.5 winner spot + 0.5 loser short perp (gross 1, net 0).
    mech 'mr': log(BTC/ETH) ratio mean reversion, z over L closes ending at
      the prior day's close. Positive state = long BTC / short ETH.
      mode 'banded': enter |z|>=2 (state = -sign(z)), exit |z|<=1, hold between.
      mode 'zfade': flat when |z|<=1, else state = -sign(z) scaled by
        min(|z|/2, 1).
    """
    w = {k: 0.0 for k in ("btc_spot", "eth_spot", "btc_perp", "eth_perp")}
    if mech == "rm":
        if r_btc is None or r_eth is None:
            return w, None
        if r_btc >= r_eth:  # BTC won the trailing window
            w["btc_spot"] += 0.5 if form == "dn" else 1.0
            if form == "dn":
                w["eth_perp"] -= 0.5
        else:
            w["eth_spot"] += 0.5 if form == "dn" else 1.0
            if form == "dn":
                w["btc_perp"] -= 0.5
        return w, None
    if mech == "mr":
        if ratio_z is None:
            return w, None
        z = ratio_z
        state = prev_state if prev_state is not None else 0.0
        if mode == "banded":
            if z <= -2.0:
                state = 1.0    # BTC weak vs ETH: expect reversion up -> long BTC
            elif z >= 2.0:
                state = -1.0
            elif abs(z) <= 1.0:
                state = 0.0
            # between bands: hold hysteresis state
        else:  # zfade
            if abs(z) <= 1.0:
                state = 0.0
            else:
                state = -1.0 if z > 0 else 1.0
        if state == 0.0:
            return w, 0.0
        scale = 1.0 if mode == "banded" else min(abs(z) / 2.0, 1.0)
        if state > 0:
            w["btc_spot"] = 0.5 * scale
            w["eth_perp"] = -0.5 * scale
        else:
            w["eth_spot"] = 0.5 * scale
            w["btc_perp"] = -0.5 * scale
        return w, state
    raise ValueError("unknown mechanism " + str(mech))


def simulate(trial, params, spot, perp, fund_btc, fund_eth, dseq, stress=1.0):
    """Run the ledger day by day. Returns daily rows.

    NAV chain: day d return = NAV(end of day d) / NAV(post-fill day d) - 1,
    where end of day d is the NAV just before the day d+1 trade (after any
    funding event settling at 00:00 d+1, charged to day d's position).
    Costs are charged inside the chain via post-fill NAV (cash paid).
    """
    ledger = Ledger(1.0)
    prev_state = None
    mech = params["mechanism"]
    rows = []
    for i, d in enumerate(dseq):
        o_prices = _open_prices(d, spot, perp)
        c_prices = _close_prices(d, spot, perp)
        if o_prices is None or c_prices is None:
            if rows:
                rows.append({"date": rows[-1]["date"] + "?gap", "ret": None,
                             "turnover": 0.0, "cost": 0.0, "net": 0.0,
                             "gross": 0.0})
            continue
        nav_open = ledger.mark(o_prices)
        # signal strictly from closes through day i-1
        tgt_w, new_state = _signal(mech, params, dseq, i, spot, prev_state)
        prev_state = new_state
        tgt_units = {s: tgt_w[s] * nav_open for s in tgt_w}
        turnover, cost = ledger.trade_to(tgt_units, o_prices, stress)
        nav_after_fill = ledger.mark(o_prices)
        # funding: (events, slot) pairs; event owns = day whose position held
        # the settled interval (00:00:xx events charged to previous day)
        for events, slot in ((fund_btc, "btc_perp"), (fund_eth, "eth_perp")):
            sym = "BTCUSDT" if slot == "btc_perp" else "ETHUSDT"
            for ts, rate in events:
                dt = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc)
                dkey = dt.date()
                if dt.hour == 0 and dt.minute == 0 and dt.second == 0:
                    dkey = dkey - timedelta(days=1)
                if dkey.isoformat() != d:
                    continue
                mp = perp_mark_lookup(perp[sym], ts)
                if mp is not None:
                    ledger.funding_event(slot, rate, mp)
        is_last = (i == len(dseq) - 1)
        if is_last:
            nav_end = ledger.mark(c_prices)
        else:
            o2 = _open_prices(dseq[i + 1], spot, perp)
            nav_end = ledger.mark(o2) if o2 is not None else ledger.mark(c_prices)
        nav_for_frac = abs(nav_after_fill)
        rows.append({
            "date": d, "ret": nav_end / nav_after_fill - 1.0,
            "turnover": turnover / nav_for_frac if nav_for_frac > 0 else 0.0,
            "cost": cost / nav_for_frac if nav_for_frac > 0 else 0.0,
            "net": _net(ledger, c_prices) / nav_for_frac if nav_for_frac > 0 else 0.0,
            "gross": _gross(ledger, c_prices) / nav_for_frac if nav_for_frac > 0 else 0.0,
        })
    return rows


OPEN_CACHE = {}


def _open_prices(d, spot, perp):
    key = ("op", d)
    if key in OPEN_CACHE:
        return OPEN_CACHE[key]
    out = {}
    for sym, slot_map in (("BTCUSDT", ("btc_spot", "btc_perp")),
                          ("ETHUSDT", ("eth_spot", "eth_perp"))):
        b = spot[sym].get(d)
        pb = perp[sym].get(d)
        if b is None or pb is None:
            return None
        out[slot_map[0]] = b[0]
        out[slot_map[1]] = pb[0]
    OPEN_CACHE[key] = out
    return out


def _close_prices(d, spot, perp):
    key = ("cl", d)
    if key in OPEN_CACHE:
        return OPEN_CACHE[key]
    out = {}
    for sym, slot_map in (("BTCUSDT", ("btc_spot", "btc_perp")),
                          ("ETHUSDT", ("eth_spot", "eth_perp"))):
        b = spot[sym].get(d)
        pb = perp[sym].get(d)
        if b is None or pb is None:
            return None
        out[slot_map[0]] = b[1]
        out[slot_map[1]] = pb[1]
    OPEN_CACHE[key] = out
    return out


def _net(ledger, prices):
    return (ledger.units["btc_spot"] * prices["btc_spot"]
            + ledger.units["eth_spot"] * prices["eth_spot"]
            + ledger.units["btc_perp"] * prices["btc_perp"]
            + ledger.units["eth_perp"] * prices["eth_perp"])




def _gross(ledger, prices):
    return sum(abs(ledger.units[s]) * prices[s] for s in ledger.units)


def _signal(mech, params, dseq, i, spot, prev_state):
    """Signal strictly from data through day i-1's close."""
    flat = ({k: 0.0 for k in ("btc_spot", "eth_spot", "btc_perp", "eth_perp")}, None)
    if i == 0:
        return flat
    L = params["lookback"]
    if L <= 0 or i < L:
        return flat
    r_btc = r_eth = ratio_z = None
    cb0, ce0 = _close_c(dseq[i - L], spot)
    cb1, ce1 = _close_c(dseq[i - 1], spot)
    if cb0 and cb1 and ce0 and ce1:
        r_btc = cb1 / cb0 - 1.0
        r_eth = ce1 / ce0 - 1.0
    if mech == "mr":
        logs = []
        for dd in dseq[i - L:i]:  # L closes ending day i-1, all <= prev close
            cb, ce = _close_c(dd, spot)
            if cb and ce and ce > 0:
                logs.append(math.log(cb / ce))
        if len(logs) == L and L >= 3:
            mu = sum(logs) / L
            var = sum((x - mu) ** 2 for x in logs) / (L - 1)
            sd = math.sqrt(var)
            if sd > 1e-12:
                ratio_z = (logs[-1] - mu) / sd
    return weights_for_rule(mech, params["form"], params["mode"],
                            r_btc, r_eth, ratio_z, prev_state)


def _close_c(d, spot):
    b = spot["BTCUSDT"].get(d)
    e = spot["ETHUSDT"].get(d)
    c_b = b[1] if b else None
    c_e = e[1] if e else None
    return c_b, c_e


# ---------------------------------------------------------------- benchmarks
def bench_ew_bah(spot, perp, dseq):
    """Equal-weight BTC/ETH spot buy-and-hold; 12 bps entry per leg once."""
    cash = 1.0
    qb = qe = 0.0
    entered = False
    rows = []
    for i, d in enumerate(dseq):
        o = _open_prices(d, spot, perp)
        c = _close_prices(d, spot, perp)
        if o is None:
            continue
        if not entered:
            cash -= 0.5 * SPOT_BPS
            cash -= 0.5 * SPOT_BPS
            qb = (cash * 0.5) / o["btc_spot"]
            qe = (cash * 0.5) / o["eth_spot"]
            cash = 0.0
            entered = True
        last = (i == len(dseq) - 1)
        px = c if last else _open_prices(dseq[i + 1], spot, perp)
        if px is None:
            px = c
        nav = cash + qb * px["btc_spot"] + qe * px["eth_spot"]
        rows.append({"date": d, "ret": None})
        if len(rows) >= 2:
            prev_nav = rows[-2]["nav"]
            rows[-1]["ret"] = nav / prev_nav - 1.0
        else:
            rows[-1]["ret"] = nav / 1.0 - 1.0
        rows[-1]["nav"] = nav
    return rows


def bench_static_dn(spot, perp, fund_eth, dseq):
    """Static dollar-neutral pair: long 0.5 BTC spot, short 0.5 ETH perp,
    entered once at first fill, never rebalanced. Funding accrues on the
    short ETH perp quantity (short receives funding when rate is negative).
    Quantities fixed at entry (deliberately unrebalanced).
    NAV chain: day d return = NAV(next open or close) / NAV(post-fill) - 1."""
    cash = 1.0
    entered = False
    qb = 0.0  # BTC spot quantity (long)
    q = 0.0   # ETH perp quantity (negative = short)
    rows = []
    for i, d in enumerate(dseq):
        o = _open_prices(d, spot, perp)
        c = _close_prices(d, spot, perp)
        if o is None:
            continue
        nav_pre = 1.0
        if not entered:
            qb = 0.5 / o["btc_spot"]
            q = -0.5 / o["eth_perp"]
            cash -= 0.5 * SPOT_BPS      # 12 bps on the spot leg notional
            cash -= 0.5 * PERP_BPS      # 7 bps on the perp leg notional
            cash -= 0.5                 # buy the spot leg; 1.0 collateral covers 0.5 short notional
            entered = True
        else:
            nav_pre = (cash + qb * o["btc_spot"] + q * o["eth_perp"]) if rows else 1.0
        for ts, rate in fund_eth:
            dt = datetime.fromtimestamp(ts / 1000.0, tz=timezone.utc)
            dkey = dt.date().isoformat()
            if dt.hour == 0 and dt.minute == 0 and dt.second == 0:
                dkey = (dt.date() - timedelta(days=1)).isoformat()
            if dkey != d:
                continue
            mp = perp_mark_lookup(perp["ETHUSDT"], ts)
            if mp is not None:
                cash += -q * rate * mp
        last = (i == len(dseq) - 1)
        px = c if last else _open_prices(dseq[i + 1], spot, perp)
        if px is None:
            px = c
        nav_end = cash + qb * px["btc_spot"] + q * px["eth_perp"]
        ret = nav_end / nav_pre - 1.0
        rows.append({"date": d, "ret": ret, "nav": nav_end})
    return rows

# ---------------------------------------------------------------- metrics
def summarize(rows, benchmark=None, bench_dn=None):
    rets = [r["ret"] for r in rows if r["ret"] is not None]
    n = len(rets)
    if n == 0:
        return {"error": "no observations"}
    mu = sum(rets) / n
    var = sum((x - mu) ** 2 for x in rets) / (n - 1) if n > 1 else 0.0
    sd = math.sqrt(var)
    sharpe = (mu / sd) * math.sqrt(ANN) if sd > 0 else 0.0
    nav = 1.0
    peak = 1.0
    maxdd = 0.0
    for x in rets:
        nav *= (1.0 + x)
        peak = max(peak, nav)
        maxdd = min(maxdd, nav / peak - 1.0)
    cum = nav - 1.0
    out = {
        "observations": n,
        "ann_return": (1.0 + cum) ** (ANN / max(n, 1)) - 1.0 if cum > -1 else -1.0,
        "ann_vol": sd * math.sqrt(ANN),
        "sharpe": sharpe,
        "cum_return": cum,
        "max_drawdown": maxdd,
        "avg_daily_turnover": sum(r["turnover"] for r in rows) / max(n, 1),
        "total_cost_pct": sum(r["cost"] for r in rows),
        "avg_net_exposure": sum(r["net"] for r in rows) / max(n, 1),
        "avg_gross_exposure": sum(r["gross"] for r in rows) / max(n, 1),
    }
    if benchmark:
        b = benchmark
        both = [(r["ret"], b[i]["ret"]) for i, r in enumerate(rows)
                if r["ret"] is not None and b[i]["ret"] is not None]
        if both:
            xs = [p[1] for p in both]
            ys = [p[0] for p in both]
            mx = sum(xs) / len(xs)
            my = sum(ys) / len(ys)
            cov = sum((x - mx) * (y - my) for x, y in both)
            vx = sum((x - mx) ** 2 for x in xs)
            beta = cov / vx if vx > 0 else 0.0
            alpha_daily = my - beta * mx
            bmu = sum(xs) / len(xs)
            bvar = sum((x - bmu) ** 2 for x in xs) / (len(xs) - 1)
            bsd = math.sqrt(bvar)
            exc = [y - x for y, x in both]
            emu = sum(exc) / len(exc)
            esd = math.sqrt(sum((e - emu) ** 2 for e in exc) / (len(exc) - 1))
            out["benchmark_sharpe"] = (bmu / bsd) * math.sqrt(ANN) if bsd > 0 else 0.0
            out["benchmark_cum"] = _prod(xs) - 1.0
            out["excess_cum"] = _prod(ys) / _prod(xs) - 1.0
            out["excess_sharpe_geo"] = 0.0
            out["market_beta_vs_bench"] = beta
            out["alpha_annualized"] = alpha_daily * ANN
            out["excess_mean_daily"] = emu
            out["excess_sharpe_arith"] = (emu / esd) * math.sqrt(ANN) if esd > 0 else 0.0
    if bench_dn:
        dn = bench_dn
        both = [(r["ret"], dn[i]["ret"]) for i, r in enumerate(rows)
                if r["ret"] is not None and i < len(dn)]
        if both:
            xs = [p[1] for p in both]
            ys = [p[0] for p in both]
            out["static_dn_excess_cum"] = _prod(ys) / _prod(xs) - 1.0
    return out


def _prod(xs):
    p = 1.0
    for x in xs:
        p *= (1.0 + x)
    return p


def summarize_bench(rows):
    rets = [r["ret"] for r in rows if r["ret"] is not None]
    n = len(rets)
    mu = sum(rets) / n if n else 0.0
    var = sum((x - mu) ** 2 for x in rets) / (n - 1) if n > 1 else 0.0
    sd = math.sqrt(var)
    nav = _prod(rets)
    peak = 1.0
    dd = 0.0
    cur = 1.0
    for x in rets:
        cur *= (1.0 + x)
        peak = max(peak, cur)
        dd = min(dd, cur / peak - 1.0)
    return {
        "observations": n,
        "sharpe": (mu / sd) * math.sqrt(ANN) if sd > 0 else 0.0,
        "cum_return": nav - 1.0,
        "max_drawdown": dd,
    }


# ---------------------------------------------------------------- trials
TRIALS = [
    {"id": "T1_rm30_cash", "mechanism": "relative_momentum", "form": "cash",
     "lookback": 30, "mode": None},
    {"id": "T2_rm90_cash", "mechanism": "relative_momentum", "form": "cash",
     "lookback": 90, "mode": None},
    {"id": "T3_rm30_dn", "mechanism": "relative_momentum", "form": "dn",
     "lookback": 30, "mode": None},
    {"id": "T4_rm90_dn", "mechanism": "relative_momentum", "form": "dn",
     "lookback": 90, "mode": None},
    {"id": "T5_mr20_dn_banded", "mechanism": "ratio_mean_reversion",
     "form": "dn", "lookback": 20, "mode": "banded"},
    {"id": "T6_mr60_dn_banded", "mechanism": "ratio_mean_reversion",
     "form": "dn", "lookback": 60, "mode": "banded"},
    {"id": "T7_mr20_dn_zfade", "mechanism": "ratio_mean_reversion",
     "form": "dn", "lookback": 20, "mode": "zfade"},
    {"id": "T8_mr60_dn_zfade", "mechanism": "ratio_mean_reversion",
     "form": "dn", "lookback": 60, "mode": "zfade"},
]


def rule_params(trial):
    mech = "rm" if trial["mechanism"] == "relative_momentum" else "mr"
    return {"mechanism": mech, "lookback": trial["lookback"],
            "form": trial["form"], "mode": trial["mode"]}


# ---------------------------------------------------------------- main
def parse_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--rule", default=None,
                    help="trial id; default = frozen selection in strategy.json")
    ap.add_argument("--stress", type=float, default=1.0,
                    help="cost multiplier (1.0 or 2.0)")
    ap.add_argument("--all-trials", action="store_true",
                    help="evaluate the whole predeclared trial table (dev window "
                         "reproduction only; refuses end > 2023-12-31)")
    ap.add_argument("--with-benchmarks", action="store_true",
                    help="also emit benchmark return series in metrics.json")
    return ap.parse_args(argv)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    start = date.fromisoformat(args.start)
    end = date.fromisoformat(args.end)
    log = []
    if end > WAVE_ONE_LAST_DAY and os.environ.get("QUANT_HOLDOUT_UNLOCK") != "1":
        sys.stderr.write(
            "refusing --end %s beyond wave-one limit 2024-12-31; set "
            "QUANT_HOLDOUT_UNLOCK=1 to unlock the embargoed window\n" % end)
        return 2
    if args.all_trials and end > DEV_END and os.environ.get("QUANT_HOLDOUT_UNLOCK") != "1":
        sys.stderr.write("--all-trials is dev-window reproduction only; "
                         "validation data may be compared once per frozen rule. Use --rule.\n")
        return 2
    if start < date(2020, 8, 1):
        sys.stderr.write("no data before 2020-08-01\n")
        return 2

    os.makedirs(args.output_dir, exist_ok=True)
    _spot_raw = {s: load_hourly(s, "spot", start, end, log) for s in ("BTCUSDT", "ETHUSDT")}
    _perp_raw = {s: load_hourly(s, "perp", start, end, log) for s in ("BTCUSDT", "ETHUSDT")}
    fund_btc = load_funding("BTCUSDT", start, end, log)
    fund_eth = load_funding("ETHUSDT", start, end, log)
    spot = {s: daily_bars(_spot_raw[s], log, "spot-" + s) for s in ("BTCUSDT", "ETHUSDT")}
    perp = {s: daily_bars(_perp_raw[s], log, "perp-" + s) for s in ("BTCUSDT", "ETHUSDT")}

    all_days = sorted(set(spot["BTCUSDT"]) & set(spot["ETHUSDT"])
                      & set(perp["BTCUSDT"]) & set(perp["ETHUSDT"]))
    dseq = [d for d in all_days if start.isoformat() <= d <= end.isoformat()]
    log.append("days=%d first=%s last=%s" % (len(dseq), dseq[0] if dseq else "-",
                                             dseq[-1] if dseq else "-"))

    results = {"log": log, "window": {"start": args.start, "end": args.end,
                                      "days": len(dseq)}}
    if args.all_trials:
        table = []
        for t in TRIALS:
            rows = simulate(t["id"], rule_params(t), spot, perp, fund_btc,
                            fund_eth, dseq, stress=args.stress)
            bah = bench_ew_bah(spot, perp, dseq)
            m = summarize(rows, benchmark=bah)
            table.append({"id": t["id"], "trial": t, "stress": args.stress,
                          "metrics": m})
        results["trials"] = table
    else:
        rule = args.rule
        if rule is None:
            rule = frozen_rule_id()
        t = next((x for x in TRIALS if x["id"] == rule), None)
        if t is None:
            sys.stderr.write("unknown rule %s\n" % rule)
            return 2
        rows = simulate(t["id"], rule_params(t), spot, perp, fund_btc,
                        fund_eth, dseq, stress=args.stress)
        bah = bench_ew_bah(spot, perp, dseq)
        dn = bench_static_dn(spot, perp, fund_eth, dseq)
        m = summarize(rows, benchmark=bah, bench_dn=dn)
        results["rule"] = t
        results["metrics"] = m
        results["benchmark_ew_bah"] = summarize_bench(bah)
        results["benchmark_static_dn"] = summarize_bench(dn)
        if args.with_benchmarks:
            results["series"] = {
                "strategy": [{"date": r["date"], "ret": r["ret"]} for r in rows],
                "ew_bah": [{"date": r["date"], "ret": r["ret"]} for r in bah],
            }
        with open(os.path.join(args.output_dir, "returns.csv"), "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["date", "strategy_return", "benchmark_return",
                         "position", "turnover", "cost"])
            for i, r in enumerate(rows):
                bret = bah[i]["ret"] if i < len(bah) and bah[i]["date"] == r["date"] else ""
                wr.writerow([r["date"],
                             "" if r["ret"] is None else "%.10f" % r["ret"],
                             "" if bret is None else "%.10f" % bret,
                             "%.6f" % r["net"], "%.6f" % r["turnover"],
                             "%.10f" % r["cost"]])
    with open(os.path.join(args.output_dir, "metrics.json"), "w") as fh:
        json.dump(results, fh, indent=1, sort_keys=True)
    return 0


FROZEN_CACHE = {}


def frozen_rule_id():
    path = os.path.join(HERE, "strategy.json")
    if path in FROZEN_CACHE:
        return FROZEN_CACHE[path]
    if not os.path.exists(path):
        raise RuntimeError("strategy.json with frozen selected_trial_id not "
                           "found; pass --rule for reproduction runs")
    with open(path) as fh:
        doc = json.load(fh)
    FROZEN_CACHE[path] = doc["selected_trial_id"]
    return doc["selected_trial_id"]


if __name__ == "__main__":
    sys.exit(main())

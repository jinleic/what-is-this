# Deribit Taker-Flow GEX — Source Verification Notes (observed 2026-08-30)

All observations below are from live, unauthenticated calls made on 2026-08-30 from this
workspace, cross-checked against the official documentation at docs.deribit.com (mintlify
mirror of Deribit API docs). Production endpoints:

- REST:  https://www.deribit.com/api/v2/...
- WS:    wss://www.deribit.com/ws/api/v2
- (History derivative host `history.deribit.com/api/v2` also serves trades unauthenticated.)

## 1. Instrument metadata (strike, expiry, kind, multiplier)

Verified live: `public/get_instruments?currency=BTC&kind=option&expired=false` → 1020 open
BTC options; per-record fields include:

- `instrument_name` e.g. `BTC-11SEP26-60000-C`; `option_type`: `"call"|"put"`;
  `strike` (BTC-USD units, e.g. 60000.0);
- `expiration_timestamp` (ms), `creation_timestamp` (ms);
- `contract_size` = 1.0 for BTC options (DK: options deliver **base currency** — the
  contract size is in BTC/ETH coins; per official docs `public/get_contract_size`:
  "value of one contract ... in USD for inverse futures/perps, **base currency coin for
  options**"); `lot_size` (integer), `tick_size`, `block_trade_tick_size`;
- `settlement_currency` = BTC (inverse options), `price_index` = `btc_usd`,
  `underlying_type` = `crypto`, `instrument_type` = `reversed`;
- `maker_commission` / `taker_commission` = 0.0003 (0.03%);
- `is_active`, `state` = `open`.

`public/get_instrument?instrument_name=BTC-11SEP26-60000-C` returns the same record
solo; verified `strike=60000.0`, `option_type=call`, `expiration_timestamp=1789113600000`,
`contract_size=1.0`, `lot_size=1`, `tick_size=0.0001`, `settlement_currency=BTC`,
`settlement_period=week`.

Docs: instrument `amount`/`contracts` units for options are the **underlying coin**
(BTC or ETH), NOT USD (unlike inverse futures/perps where amount is USD).

## 2. Index price

Verified live: `public/get_index_price?index_name=btc_usd` →
`{"index_price": 78153.81, "estimated_delivery_price": 78153.81}` (USD).
Docs list `deribit_price_index.{index_name}` WS channel for pushes; the option ticker
pushes `index_price` and `underlying_price` fields directly.

## 3. Trade semantics (taker direction)

Verified live on `public/get_last_trades_by_currency?currency=BTC&kind=option` — fields:
`timestamp` (ms), `instrument_name`, `direction` (`buy`/`sell`), `amount` (coin),
`contracts` (number of contracts), `price` (option premium in **BTC**, the quote
currency of the inverse option), `mark_price`, `iv` (%), `index_price` (USD at trade
time), `trade_id` (string, unique per currency), `trade_seq` (per instrument),
`tick_direction` (0=Plus,1=Zero-Plus,2=Minus,3=Zero-Minus).

Official docs (`trades.(kind).(currency).(interval)` channel, AsyncAPI schema) mark
`direction` as the **taker** direction: `buy` = taker bought the option (aggressive
buyer, dealer sold), `sell` = taker sold the option (aggressive seller, dealer bought).
Optional fields: `liquidation` (`M`/`T`/`MT`), `block_trade_id`,
`block_trade_leg_count`, `combo_id`, `combo_trade_id`, `starbase_match_id`,
`block_rfq_id`, `starbase_timestamp` (ns).

**Taker sign convention for dealer-mirror GEX** (recorded, preregistered):
signed coin amount = +amount if direction=buy, −amount if sell.
GEX contribution per trade = signed_amount(coin) × gamma × S² × 0.01 / K...
see run.py: recorded exactly as `signed_gamma_usd` with:
    signed_usd = sign × amount Coins × price(BTC) × ... — units reconciled in §4.

## 4. Greeks and units

Verified live: `public/ticker?instrument_name=BTC-11SEP26-60000-C` →
`greeks = {"delta": 0.99245, "gamma": 0.0, "vega": 2.93958, "theta": -7.68748, "rho": 19.3006}`.

Docs (`ticker.{instrument}.{interval}` channel / `public/ticker`): greeks object is
options-only, "calculated using standard Black Scholes without adjustments" with these
scale caveats stated by Deribit:

- `gamma`: **per 1 USD change in underlying** (standard Black-Scholes ∂Δ/∂S, S in USD).
- `vega` is per 1 USD absolute underlying move (not per vol pt); verified scale:
  nearer-the-money BTC option ticker gave `vega: 0.00038`-style values in BTC units
  (premium-coin units). DK: vega is in premium units per USD move — check live value
  in the smoke run; not load-bearing for GEX (we use gamma only).
- `theta`: **minimum of (1-day theta, lifetime theta)** — nonstandard near expiry;
  recentred in capture because it is NOT used by GEX.
- `delta`: Black–Scholes delta (not Net Transaction Delta; DeltaTotal = NTD + futures
  BTC position applies only to account-summary objects).

**Gamma units** (preregistered reconciliation for inverse coin-margined options):
standard Black-Scholes gamma in 1/USD², i.e. ∂ΔBTC/∂S_USD per $1 move... per 1 USD:
`gamma_sched = gamma` where `Δ` measured in **BTC delta per 1 USD underlying change**.
Dealer-mirror USD gamma per contract:
    GEX_usd_per_contract = gamma × S²  [USD of delta per 1 USD move], with S=index
    price in USD. With `amount` in coins (contracts × contract_size = same number since
    contract_size=1), total taker-signed gamma in USD per 1% move =
    signed_amount(coins) × gamma × S² × 0.01.
This matches the standard convention
    GEX$ = Σ_signed contracts × contract_multiplier × gamma × S² × 0.01
with contract_multiplier = 1 coin (contract_size) for Deribit options, and amounts in
coins. Positions in open_interest (coin units) at `public/get_book_summary_by_currency`
also use coin units per docs.

## 5. WS channel semantics (public, unauthenticated)

Observed live (stdlib WS client):

- `wss://www.deribit.com/ws/api/v2` — handshake `HTTP/1.1 101 Switching Protocols` OK.
- `public/subscribe {"channels":["trades.option-btc_usd","ticker.X.raw"]}` → **error
  13778 `raw_subscriptions_not_available_for_unauthorized`** whenever ANY channel in the
  batch uses `.raw` (whole batch fails: ack returns `error`, no channels subscribed).
- Non-raw batch subscribe with DASH-form names
  `["trades.future-btc_usd.100ms","trades.option-btc_usd.100ms"]` → acked
  **`result: []` — a REJECTION in practice**: the server acks without error but
  subscribes nothing (no pushes ever arrive on those names). Canonical dot-form
  names echo the full requested set in `result`; readiness requires that
  ack-set equality, and `result: []` is never success.
- High-frequency push proof: `ticker.BTC-PERPETUAL.100ms` pushes ~3–5 messages/sec with
  full ticker payload (parse verified). No option trade occurred during the first
  60–90 s windows; live option trades occur in bursts (REST shows several/minute).
- Docs: `trades.{kind}.{currency}.{interval}` intervals `raw|100ms|agg2`; raw
  unauthorized. Channel names in push `params.channel`.
- Docs: `public/set_heartbeat` enables `heartbeat` notifications → client must respond
  `public/test`; serves liveness. Server also sends RFC6455 WS-level pings (opcode 9)
  which must be ponged (opcode 10).
- Instrument metadata refresh: docs advise seeding from `public/get_instruments` and
  tracking `instrument.creation.{kind}.{currency}` / `instrument.state.{kind}.{currency}`
  instead of polling (get_instruments rate limit is much lower than other endpoints).
- Chain greeks alternative: `markprice.options.{index_name}` streams mark price + IV for
  the chain but **NOT gamma**; gamma only comes from `ticker.{instrument}.{interval}`
  (or `public/ticker` REST) per instrument, or `incremental_ticker` (≤1 Hz/instrument).

## 6. Fields that CANNOT be reconstructed from public streams

- `.raw` (1 ms) aggregation for any channel without authorization; public minimum is
  100ms aggregation (or agg2). GEX accounting tolerates this (no regime claim before
  four registered weeks anyway), but 1 ms interleaving is unavailable.
- `user.*` channels (own orders/trades/etc.) — out of scope (research-only, no account).
- Historical trade backfill: `public/get_last_trades_by_currency_and_time` (count ≤1000,
  paging by advancing timestamp; sequence paging `start_seq/end_seq` on
  `get_last_trades_by_instrument`) — verified semantic from docs; depends on venue
  retention limits (docs warn of limited lookback; not separable from full history).
- Starbase feed has NO greeks/IV (docs) — only standard ticker channels carry gamma.

## 7. Trailing facts

- `public/get_currencies` works unauthenticated; BTC and ETH option markets exist
  (1020 BTC, 862 ETH open options at observation time; expiries daily + weekly + far).
- `public/get_chart_trades` does NOT exist (method not found); candle endpoint is
  `public/get_tradingview_chart_data`.
- WS subscribe semantics (CORRECTED in section 8): an ack carrying `result: []`
  means **NOTHING was subscribed** — for the dash-form channel name the server
  acks without error yet the returned set is empty. It is therefore a
  REJECTION, never treated as success; readiness requires the ack result to be
  exactly the requested channel set. Error code 13778 = raw-only rejection.
- Book summary by instrument includes `mark_iv` (%) 62.04, `underlying_index` =
  `BTC-11SEP26` (underlying expiry), `interest_rate`, OI in coins.


## 8. Corrective round (parent-supervised, 2026-08-30)

Additional live-verified facts that changed the collector configuration:

- **Canonical trade channel syntax**: `trades.option.BTC.100ms` /
  `trades.option.ETH.100ms` acked with `result` echoing exactly those channel
  names. (The dash form `trades.option-btc_usd.100ms` also acks but returns an
  empty `result` array, i.e. it silently subscribes nothing — which is why
  readiness now requires ack-set equality.)
- **`public/set_heartbeat` REQUIRES `interval`**: sending `params: {}` returns
  `-32602 Invalid params {"reason":"value required","param":"interval"}`.
  `params: {"interval": 10}` returns `result: "ok"`. `public/disable_heartbeat`
  returns `result: "ok"`. Client answers `heartbeat` notifications with
  `public/test`.
- **Sustainable default**: only 4 channels (2 option-trade + 2 index). The
  chain-wide `ticker.{instrument}.agg2` fan-out (1,882 channels) wrote ~166 MiB
  in 60 s and would exhaust a 256 MiB/day cap within minutes; it is now behind
  `--capture-tickers` and is NOT default. Index pushes are ~1–2/s per currency.
- **Live parse proof on the sustainable default**: 30 real option trades in a
  75 s window with `direction` (taker), `amount` (coin), `price` (coin premium),
  `index_price` (USD), `iv` (%), `trade_id`, `trade_seq`; 141 index records
  (`btc_usd`, `eth_usd`); zero server errors after the heartbeat fix.
- **Gamma without ticker capture**: replay computes standard Black-Scholes gamma
  from the trade's own embedded `iv` (%), plus `strike` and
  `expiration_timestamp` from instrument metadata and S = trade `index_price`,
  labeled `greeks_source = "bs_from_trade_iv"` (model-derived, NOT a venue
  value). Venue gamma from captured tickers, when present, is labeled
  `greeks_source = "deribit_ticker"`. Missing `iv` or missing strike/expiry
  fails closed: the trade is excluded and counted in coverage.
- **Byte accounting**: exact encoded payload bytes, `flush` + `os.fsync` per
  record, startup reconciliation over every `dayYYYYMMDD/` directory, active-day
  bytes loaded at first rotate, and the rotation-open event itself checked
  against both caps so a restart cannot bypass them.


## 9. Review-driven repairs (2026-08-30, post independent review)

Review verdict: LAUNCH BLOCKED. Repairs implemented (all no-network proof-checked
in `proofs_capture.py`, 20/20 pass):

- Rotation: Sink.write rotates BEFORE every write whenever the capture UTC day
  changes (not only at first use); the day_rotation_open record itself is
  checked against both caps so it cannot bypass them.
- Day-scoped metadata: on day rollover `_persisted_instruments` resets, so the
  new day's instruments.jsonl always carries strike/expiry/contract_size
  (needed for the default model-derived BS gamma). Same-day reconnects never
  rewrite unchanged rows (rehydrated ids + `_instrument_persist_day`).
- Exact ack equality: subscribe batches succeed ONLY when the ack result is
  set-equal to the request. Empty, partial, superset, errored, and timed-out
  acks all fail and force a reconnect (never a live but unready loop).
- Readiness: DERIBIT_CAPTURE_READY prints ONLY after every batch ack is exact
  AND set_heartbeat (interval required, result=ok) is acked; the receive loop
  is entered only then. Any failure reconnects with capped backoff.
- Backoff resets only after a stable ready session (not after TCP setup).
- All waits (reconnect sleep, REST retries, ack waits) are SIGTERM-bounded via
  stop_aware_sleep / request_stop; observed SIGTERM latency 0 ms.
- Dossier wording: raw JSONL is AT-LEAST-ONCE across restarts with replay dedup
  authoritative (run.py deterministic ids). All stale dash-channel/result:[]/
  ticker-default/lifecycle claims replaced with the canonical 4-channel
  configuration and lifecycle semantics.

Proof scenarios (deterministic, no network):
1  midnight rotate-before-write; old day untouched; new day gets rotation record + write
2  day cap fail-closed
3  cumulative cap across restart; no debris
4  rotation record itself blocked by pre-existing day bytes
5  ack matrix: exact / empty / partial / superset / error
6  heartbeat setup error blocks readiness
7  failed subscribe -> no readiness; backoff escalates 1->2->4
8  SIGTERM during backoff sleep: wakes in 0 ms
9  SIGTERM pending REST: immediate abort with 'stop requested'
10 day-scoped metadata: same-day reconnect writes nothing; rollover re-persists

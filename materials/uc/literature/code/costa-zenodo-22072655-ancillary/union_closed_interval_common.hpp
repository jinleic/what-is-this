#pragma once
#if __has_include(<gmp.h>)
#  include <gmp.h>
#else
typedef unsigned long int mp_limb_t;
#endif
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::radix == 2 &&
              std::numeric_limits<double>::digits == 53 &&
              std::numeric_limits<double>::min_exponent == -1021 &&
              std::numeric_limits<double>::max_exponent == 1024 &&
              sizeof(double) == 8,
              "This verifier requires IEEE-754 binary64 double semantics.");

#if __has_include(<mpfr.h>)
#  include <mpfr.h>
#else
extern "C" {
typedef long int mpfr_prec_t;
typedef long int mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum { MPFR_RNDN=0, MPFR_RNDZ=1, MPFR_RNDU=2, MPFR_RNDD=3,
               MPFR_RNDA=4, MPFR_RNDF=5, MPFR_RNDNA=-1 } mpfr_rnd_t;
typedef struct { mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign;
                 mpfr_exp_t _mpfr_exp; mp_limb_t *_mpfr_d; } __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct *mpfr_ptr;
typedef const __mpfr_struct *mpfr_srcptr;
void mpfr_init2(mpfr_ptr, mpfr_prec_t); void mpfr_clear(mpfr_ptr);
int mpfr_set_d(mpfr_ptr,double,mpfr_rnd_t);
int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t);
const char* mpfr_get_version(void);
}
#endif

double down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
double up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
struct I {
    double lo,hi;
    I():lo(0),hi(0){}
    I(double x):lo(x),hi(x){}
    I(double a,double b):lo(a),hi(b){if(a>b)throw std::runtime_error("bad interval");}
    double lower()const{return lo;} double upper()const{return hi;}
};
I decimal_point(double x){return I(down(x),up(x));}
I decimal_interval(double a,double b){return I(down(a),up(b));}
double report_lower(double x){return down(x);}
double report_upper(double x){return up(x);}
std::vector<double> partition(double a,double b,int n){
    std::vector<double> z(n+1); z.front()=a; z.back()=b;
    for(int i=1;i<n;++i) z[i]=a+(b-a)*i/n;
    return z;
}
I operator+(const I&a,const I&b){return I(down(a.lo+b.lo),up(a.hi+b.hi));}
I operator-(const I&a,const I&b){return I(down(a.lo-b.hi),up(a.hi-b.lo));}
I operator-(const I&a){return I(down(-a.hi),up(-a.lo));}
I operator*(const I&a,const I&b){
    double p[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    return I(down(*std::min_element(p,p+4)),up(*std::max_element(p,p+4)));
}
I reciprocal(const I&b){
    if(b.lo<=0&&b.hi>=0)throw std::runtime_error("division by zero interval");
    double p=1.0/b.lo,q=1.0/b.hi;
    return I(down(std::min(p,q)),up(std::max(p,q)));
}
I operator/(const I&a,const I&b){return a*reciprocal(b);}
static constexpr mpfr_prec_t MPFR_PREC = 192;

enum class TF { LOG, EXP, SQRT };
struct MPFRWorkspace {
    mpfr_t a,b;
    MPFRWorkspace(){ mpfr_init2(a,MPFR_PREC); mpfr_init2(b,MPFR_PREC); }
    ~MPFRWorkspace(){ mpfr_clear(a); mpfr_clear(b); }
};
double mpfr_eval(double x, TF f, mpfr_rnd_t rnd) {
    static thread_local MPFRWorkspace w;
    mpfr_set_d(w.a,x,MPFR_RNDN);
    if(f==TF::LOG) mpfr_log(w.b,w.a,rnd);
    else if(f==TF::EXP) mpfr_exp(w.b,w.a,rnd);
    else mpfr_sqrt(w.b,w.a,rnd);
    return mpfr_get_d(w.b,rnd);
}
double a4_cut() {
    // sqrt(1/2)=1/sqrt(2), with the final conversion rounded upward.
    // Thus the binary64 cut lies on the zero-branch side of the exact cut.
    static const double q=mpfr_eval(0.5,TF::SQRT,MPFR_RNDU);
    return q;
}
I ilog(const I& x) {
    if(x.lower()<=0) throw std::runtime_error("log nonpositive");
    return I(mpfr_eval(x.lower(),TF::LOG,MPFR_RNDD),
             mpfr_eval(x.upper(),TF::LOG,MPFR_RNDU));
}
I iexp(const I& x) {
    return I(mpfr_eval(x.lower(),TF::EXP,MPFR_RNDD),
             mpfr_eval(x.upper(),TF::EXP,MPFR_RNDU));
}
I isqrt(const I& x) {
    double lo=x.lower();
    if(lo<0 && x.upper()>=0) lo=0;
    if(lo<0) throw std::runtime_error("sqrt negative");
    return I(mpfr_eval(lo,TF::SQRT,MPFR_RNDD),
             mpfr_eval(x.upper(),TF::SQRT,MPFR_RNDU));
}
I sqr(const I& x) {
    if(x.lower()>=0) return x*x;
    if(x.upper()<=0) return (-x)*(-x);
    const double hi=std::max(x.lower()*x.lower(),x.upper()*x.upper());
    return I(0.0,std::nextafter(hi,std::numeric_limits<double>::infinity()));
}
I iabs(const I& x) {
    if(x.lower()>=0) return x;
    if(x.upper()<=0) return -x;
    return I(0.0,std::max(-x.lower(),x.upper()));
}

const I XSTAR=decimal_interval(0.69078759392498,0.69078759392500);
const I MSTAR=decimal_interval(0.61729091208120,0.61729091208133);
const I BETA=decimal_interval(0.1000525598620,0.1000525598640);
const I ALPHA=I(1.0)-BETA;
const I THETA=I(111.0)/I(2500.0);
const I EPSILON_CERT=I(41.0)/I(2000000.0);
const I MTARGET=MSTAR-EPSILON_CERT;
const I SIGMA=I(99999.0)/I(100000.0);
const I KAPPA=I(7.0)/I(10.0);
const I HALF(0.5);
const I LOG2(mpfr_eval(2.0,TF::LOG,MPFR_RNDD),
             mpfr_eval(2.0,TF::LOG,MPFR_RNDU));

I Hraw(const I& z){ return -z*ilog(z)-(I(1.0)-z)*ilog(I(1.0)-z); }
I Hpoint(double z){ return Hraw(I(z)); }
I Hmono(const I& z) {
    if(z.lower()<=0 || z.upper()>=1) {
        throw std::runtime_error("H outside open interval ["+std::to_string(z.lower())+","+std::to_string(z.upper())+"]");
    }
    if(z.upper()<=0.5) {
        I a=Hpoint(z.lower()), b=Hpoint(z.upper());
        return I(a.lower(),b.upper());
    }
    if(z.lower()>=0.5) {
        I a=Hpoint(z.lower()), b=Hpoint(z.upper());
        return I(b.lower(),a.upper());
    }
    I a=Hpoint(z.lower()), b=Hpoint(z.upper());
    return I(std::min(a.lower(),b.lower()),LOG2.upper());
}
I Hp(const I& z){ return ilog((I(1.0)-z)/z); }

I gxraw(const I&x){return x*x*(I(1.0)+(I(1.0)-x)*(I(1.0)-x));}
I gx(const I&x){
    I a=gxraw(I(x.lower())), b=gxraw(I(x.upper()));
    return I(a.lower(),b.upper());
}
I gpx(const I&x){return I(4.0)*x-I(6.0)*x*x+I(4.0)*x*x*x;}
I gppx(const I&x){return I(4.0)-I(12.0)*x+I(12.0)*x*x;}
I Hx(const I&x){return Hmono(x);}

// Example 4 in complemented coordinates.  These functions are called only
// on intervals contained in one of the three analytic branches.
int a4_branch(const I& x) {
    const double q=a4_cut();
    if(x.upper()<=0.5) return 0;
    if(x.lower()>=0.5 && x.upper()<=q) return 1;
    if(x.lower()>=q) return 2;
    return -1;
}
I a4(const I&x) {
    int b=a4_branch(x);
    if(b==0) return I(1.0);
    if(b==2) return I(0.0);
    if(b!=1) throw std::runtime_error("a4 branch crossing");
    I N=I(1.0)-I(2.0)*x*x;
    // The directed binary64 cut lies above the exact value 1/sqrt(2).
    // Intersect the tiny negative roundoff tail with the true range N>=0.
    if(N.lower()<0) N=I(0.0,N.upper());
    return isqrt(N/(I(2.0)*x*(I(1.0)-x)));
}
I a4p(const I&x) {
    int b=a4_branch(x);
    if(b==0 || b==2) return I(0.0);
    I A=a4(x), N=I(1.0)-I(2.0)*x*x, D=I(2.0)*x*(I(1.0)-x);
    if(N.lower()<=0) throw std::runtime_error("singular a4 derivative at branch point");
    return A*HALF*((-I(4.0)*x)/N-(I(2.0)-I(4.0)*x)/D);
}

enum class PhiBranch { ZERO_LOW, RAMP_LOW, LOW, MID, HIGH, RAMP_HIGH, ZERO_HIGH, CROSS };
PhiBranch phi_branch(const I&x) {
    const double q=181.0/256.0;
    if(x.upper()<=1.0/64.0) return PhiBranch::ZERO_LOW;
    if(x.lower()>=1.0/64.0 && x.upper()<=1.0/16.0) return PhiBranch::RAMP_LOW;
    if(x.lower()>=1.0/16.0 && x.upper()<=0.5) return PhiBranch::LOW;
    if(x.lower()>=0.5 && x.upper()<=q) return PhiBranch::MID;
    if(x.lower()>=q && x.upper()<=15.0/16.0) return PhiBranch::HIGH;
    if(x.lower()>=15.0/16.0 && x.upper()<=63.0/64.0) return PhiBranch::RAMP_HIGH;
    if(x.lower()>=63.0/64.0) return PhiBranch::ZERO_HIGH;
    return PhiBranch::CROSS;
}
I cutoff(const I&x, PhiBranch b) {
    if(b==PhiBranch::ZERO_LOW || b==PhiBranch::ZERO_HIGH) return I(0.0);
    if(b==PhiBranch::RAMP_LOW) return (x-I(1.0/64.0))/I(3.0/64.0);
    if(b==PhiBranch::RAMP_HIGH) return (I(63.0/64.0)-x)/I(3.0/64.0);
    return I(1.0);
}
I cutoffp(PhiBranch b) {
    if(b==PhiBranch::RAMP_LOW) return I(64.0)/I(3.0);
    if(b==PhiBranch::RAMP_HIGH) return -I(64.0)/I(3.0);
    return I(0.0);
}
I diag_entropy(const I&x, PhiBranch b) {
    if(b==PhiBranch::LOW || b==PhiBranch::RAMP_LOW) return Hx(x);
    if(b==PhiBranch::MID) return LOG2;
    if(b==PhiBranch::HIGH || b==PhiBranch::RAMP_HIGH)
        return LOG2*Hmono(x*x)/Hmono(I(181.0/256.0)*I(181.0/256.0));
    throw std::runtime_error("diag entropy on zero/cross branch");
}
I diag_entropyp(const I&x, PhiBranch b) {
    if(b==PhiBranch::LOW || b==PhiBranch::RAMP_LOW) return Hp(x);
    if(b==PhiBranch::MID) return I(0.0);
    if(b==PhiBranch::HIGH || b==PhiBranch::RAMP_HIGH)
        return LOG2*I(2.0)*x*Hp(x*x)/Hmono(I(181.0/256.0)*I(181.0/256.0));
    throw std::runtime_error("diag entropy derivative on zero/cross branch");
}
struct PhiData { I v, px; bool ok; };
PhiData phi_all_one_branch(const I&x) {
    PhiBranch b=phi_branch(x);
    if(b==PhiBranch::CROSS) return {I(0.0),I(0.0),false};
    if(b==PhiBranch::ZERO_LOW || b==PhiBranch::ZERO_HIGH)
        return {I(0.0),I(0.0),true};
    I HS=Hmono(XSTAR), ell=HALF*Hp(XSTAR)/HS;
    I D=diag_entropy(x,b), Dp=diag_entropyp(x,b), HX=Hx(x), HXp=Hp(x);
    I dx=x-XSTAR;
    I E=iexp(-ell*dx-KAPPA*sqr(dx));
    I base=SIGMA*isqrt(D*HX/HS)*E;
    I logder=HALF*Dp/D+HALF*HXp/HX-ell-I(2.0)*KAPPA*dx;
    I ch=cutoff(x,b), chp=cutoffp(b);
    return {base*ch,base*logder*ch+base*chp,true};
}
PhiData phi_all(const I&x) {
    if(phi_branch(x)!=PhiBranch::CROSS) return phi_all_one_branch(x);
    const double cuts[]={1.0/64.0,1.0/16.0,0.5,181.0/256.0,15.0/16.0,63.0/64.0};
    std::vector<double> p={x.lower()};
    for(double c:cuts) if(x.lower()<c && c<x.upper()) p.push_back(c);
    p.push_back(x.upper());
    bool first=true; I vh,dh;
    for(size_t i=0;i+1<p.size();++i) {
        PhiData z=phi_all_one_branch(I(p[i],p[i+1]));
        if(!z.ok) return {I(0.0),I(0.0),false};
        if(first){vh=z.v;dh=z.px;first=false;}
        else {
            vh=I(std::min(vh.lower(),z.v.lower()),std::max(vh.upper(),z.v.upper()));
            dh=I(std::min(dh.lower(),z.px.lower()),std::max(dh.upper(),z.px.upper()));
        }
    }
    return {vh,dh,true};
}

struct XYBox { double xa,xb,ya,yb; int dep; bool diag; };
bool minorant_mean_value(const XYBox& b, double& cert) {
    I X(b.xa,b.xb),Y(b.ya,b.yb);
    if(phi_branch(X)==PhiBranch::CROSS || phi_branch(Y)==PhiBranch::CROSS) return false;
    const double xc=0.5*(b.xa+b.xb), yc=0.5*(b.ya+b.yb);
    // The rounded midpoint need not be exactly equidistant from the two
    // endpoints.  Use the larger actual binary64 distance and round it up.
    const double hx=std::nextafter(std::max(xc-b.xa,b.xb-xc),INFINITY);
    const double hy=std::nextafter(std::max(yc-b.ya,b.yb-yc),INFINITY);
    auto px=phi_all(X),py=phi_all(Y),pcx=phi_all(I(xc)),pcy=phi_all(I(yc));
    I AX=a4(X),AY=a4(Y);
    I c=Y+AX*AY*(I(1.0)-Y), q=X*c;
    if(q.lower()<=0 || q.upper()>=1) return false;
    I direct=Hmono(q)-px.v*py.v;
    if(direct.lower()>0) {cert=direct.lower();return true;}
    I AXp,AYp;
    try { AXp=a4p(X); AYp=a4p(Y); }
    catch(const std::runtime_error&) { return false; }
    I qx=c+X*AXp*AY*(I(1.0)-Y);
    I qy=X*(I(1.0)+AX*(AYp*(I(1.0)-Y)-AY));
    I fX=Hp(q)*qx-px.px*py.v;
    I fY=Hp(q)*qy-px.v*py.px;
    I AXc=a4(I(xc)),AYc=a4(I(yc));
    I qc=I(xc)*(I(yc)+AXc*AYc*(I(1.0)-I(yc)));
    I fc=Hmono(qc)-pcx.v*pcy.v;
    I lower=fc-iabs(fX)*I(hx)-iabs(fY)*I(hy);
    cert=lower.lower();
    return cert>0;
}

bool certify_minorant() {
    const double q=a4_cut();
    const std::vector<double> cuts={1.0/64.0,1.0/16.0,0.5,181.0/256.0,q,15.0/16.0,63.0/64.0};
    std::vector<XYBox> stack;
    for(size_t i=0;i+1<cuts.size();++i)
        for(size_t j=i;j+1<cuts.size();++j)
            stack.push_back({cuts[i],cuts[i+1],cuts[j],cuts[j+1],0,i==j});
    long long leaves=0,splits=0; int maxd=0; double mincert=INFINITY;
    std::vector<XYBox> bad;
    while(!stack.empty()) {
        XYBox b=stack.back(); stack.pop_back(); maxd=std::max(maxd,b.dep);
        double cert=0;
        bool passed=false;
        try { passed=minorant_mean_value(b,cert); }
        catch(const std::runtime_error& e) {
            std::cerr<<"minorant exception on ["<<b.xa<<","<<b.xb<<"] x ["
                     <<b.ya<<","<<b.yb<<"]: "<<e.what()<<"\n";
            throw;
        }
        if(passed) { mincert=std::min(mincert,cert); ++leaves; continue; }
        if(b.dep>=34) { bad.push_back(b); if(bad.size()>=20) break; continue; }
        if(b.diag) {
            double m=0.5*(b.xa+b.xb);
            stack.push_back({b.xa,m,b.xa,m,b.dep+1,true});
            stack.push_back({b.xa,m,m,b.xb,b.dep+1,false});
            stack.push_back({m,b.xb,m,b.xb,b.dep+1,true});
        } else if((b.xb-b.xa)>=(b.yb-b.ya)) {
            double m=0.5*(b.xa+b.xb);
            stack.push_back({b.xa,m,b.ya,b.yb,b.dep+1,false});
            stack.push_back({m,b.xb,b.ya,b.yb,b.dep+1,false});
        } else {
            double m=0.5*(b.ya+b.yb);
            stack.push_back({b.xa,b.xb,b.ya,m,b.dep+1,false});
            stack.push_back({b.xa,b.xb,m,b.yb,b.dep+1,false});
        }
        ++splits;
    }
    std::cout<<"Example 4 product lower bound: leaves="<<leaves
             <<", splits="<<splits<<", maxdepth="<<maxd
             <<", unresolved="<<bad.size()<<", min certified lower="<<report_lower(mincert)<<"\n";
    if(!bad.empty()) {
        for(const auto& b:bad)
            std::cout<<"  bad ["<<b.xa<<","<<b.xb<<"] x ["<<b.ya<<","<<b.yb<<"]\n";
    }
    return bad.empty();
}

struct Triple { I v,p,pp; };
I G(const I&u){I x=iexp(-u);return Hmono(x)/x;}
Triple Gall(const I&u){
    I x=iexp(-u),hh=Hmono(x),hp=Hp(x);
    I v=hh/x,p=v-hp,pp=p-I(1.0)/(I(1.0)-x);
    return{v,p,pp};
}
Triple Call(const I&u){
    I x=iexp(-u),z=gx(x),J=Hmono(z),hp=Hp(z),hpp=-I(1.0)/(z*(I(1.0)-z));
    I z1=-x*gpx(x),z2=x*gpx(x)+x*x*gppx(x);
    I j1=hp*z1,j2=hpp*z1*z1+hp*z2;
    I f1=I(1.0)+HALF*j1/J,f2=HALF*(j2/J-(j1/J)*(j1/J));
    I c=isqrt(J)/x;return{c,c*f1,c*(f2+f1*f1)};
}
struct PData { I v,p; bool ok; };
PData Pall(const I&u) {
    I x=iexp(-u); auto ph=phi_all(x); if(!ph.ok) return {I(0.0),I(0.0),false};
    I p=ph.v/x;
    return {p,p-ph.px,true};
}

// On the small local rectangle used below, exp(-u) stays in the MID
// branch of phi: chi=1 and D=log(2).  There phi is C^2, so the transformed
// function P(u)=phi(exp(-u))/exp(-u) can be differentiated explicitly.
struct PData2 { I v,p,pp; bool ok; };
PData2 Pall2_mid(const I&u) {
    I x=iexp(-u);
    if(phi_branch(x)!=PhiBranch::MID) return {I(0.0),I(0.0),I(0.0),false};
    auto ph=phi_all_one_branch(x);
    I hx=Hmono(x), hp=Hp(x), hpp=-I(1.0)/(x*(I(1.0)-x));
    I ell=HALF*Hp(XSTAR)/Hmono(XSTAR), dx=x-XSTAR;
    I lx=HALF*hp/hx-ell-I(2.0)*KAPPA*dx;
    I lxp=HALF*(hpp/hx-(hp/hx)*(hp/hx))-I(2.0)*KAPPA;
    I phix=ph.v*lx, phixx=ph.v*(lx*lx+lxp);
    I p=ph.v/x, pu=p-phix, puu=pu+x*phixx;
    return {p,pu,puu,true};
}

const I LOCAL_T_LO_DEC=decimal_point(0.3520),LOCAL_T_HI_DEC=decimal_point(0.4500),LOCAL_D_HI_DEC=decimal_point(0.0050);
const double LOCAL_T_LO_OUT=LOCAL_T_LO_DEC.lower(),LOCAL_T_HI_OUT=LOCAL_T_HI_DEC.upper(),LOCAL_D_HI_OUT=LOCAL_D_HI_DEC.upper();
const double LOCAL_T_LO_IN=LOCAL_T_LO_DEC.upper(),LOCAL_T_HI_IN=LOCAL_T_HI_DEC.lower(),LOCAL_D_HI_IN=LOCAL_D_HI_DEC.lower();

struct NHessian { I v,tt,dd,td; bool ok; };
NHessian Nall_local(const I&T,const I&D) {
    I U=T+D,V=T-D;
    auto g2=Gall(I(2.0)*T),gu=Gall(U),gv=Gall(V);
    auto cu=Call(U),cv=Call(V);
    auto pu=Pall2_mid(U),pv=Pall2_mid(V);
    if(!pu.ok||!pv.ok) return {I(0.0),I(0.0),I(0.0),I(0.0),false};
    I m0=MTARGET*(I(1.0)-THETA),m4=MTARGET*THETA;
    I value=m0*(ALPHA*g2.v+BETA*cu.v*cv.v)+m4*pu.v*pv.v-HALF*(gu.v+gv.v);
    I tt=I(4.0)*m0*ALPHA*g2.pp
        +m0*BETA*(cu.pp*cv.v+I(2.0)*cu.p*cv.p+cu.v*cv.pp)
        +m4*(pu.pp*pv.v+I(2.0)*pu.p*pv.p+pu.v*pv.pp)
        -HALF*(gu.pp+gv.pp);
    I dd=m0*BETA*(cu.pp*cv.v-I(2.0)*cu.p*cv.p+cu.v*cv.pp)
        +m4*(pu.pp*pv.v-I(2.0)*pu.p*pv.p+pu.v*pv.pp)
        -HALF*(gu.pp+gv.pp);
    I td=m0*BETA*(cu.pp*cv.v-cu.v*cv.pp)
        +m4*(pu.pp*pv.v-pu.v*pv.pp)-HALF*(gu.pp-gv.pp);
    return {value,tt,dd,td,true};
}

I Nt_diagonal(const I&T) {
    auto g2=Gall(I(2.0)*T),g1=Gall(T),c=Call(T);
    auto p=Pall2_mid(T);
    if(!p.ok) throw std::runtime_error("P outside MID branch in local root check");
    I m0=MTARGET*(I(1.0)-THETA),m4=MTARGET*THETA;
    return I(2.0)*m0*(ALPHA*g2.p+BETA*c.p*c.v)
        +I(2.0)*m4*p.p*p.v-g1.p;
}

bool certify_local_hessian_and_stationary_point() {
    const int NT=40,ND=20;
    double min_tt=INFINITY,min_dd=INFINITY,min_det=INFINITY,max_trace=0.0;
    bool ok=true;
    auto tg=partition(LOCAL_T_LO_OUT,LOCAL_T_HI_OUT,NT),dg=partition(0.0,LOCAL_D_HI_OUT,ND);
    for(int i=0;i<NT;++i) for(int j=0;j<ND;++j) {
        auto z=Nall_local(I(tg[i],tg[i+1]),I(dg[j],dg[j+1]));
        if(!z.ok) return false;
        double lt=z.tt.lower(),ld=z.dd.lower();
        double mt=std::max(std::abs(z.td.lower()),std::abs(z.td.upper()));
        I det=I(lt)*I(ld)-I(mt)*I(mt);
        min_tt=std::min(min_tt,lt); min_dd=std::min(min_dd,ld);
        min_det=std::min(min_det,det.lower());
        max_trace=std::max(max_trace,z.tt.upper()+z.dd.upper());
        if(lt<=0||ld<=0||det.lower()<=0) ok=false;
    }
    const I tl=decimal_point(0.35615880),tr=decimal_point(0.35615895);
    I dl=Nt_diagonal(tl),dr=Nt_diagonal(tr);
    auto zr=Nall_local(I(tl.lower(),tr.upper()),I(0.0));
    auto zstar=Nall_local(-ilog(XSTAR),I(0.0));
    double root_value=zr.v.lower();
    ok=ok&&dl.upper()<0&&dr.lower()>0&&root_value>0&&zstar.v.lower()>0;
    std::cout<<"local rectangle: min N_tt="<<report_lower(min_tt)
             <<", min N_dd="<<report_lower(min_dd)<<", min determinant="<<report_lower(min_det)
             <<", max trace="<<report_upper(max_trace)<<"\n"
             <<"local diagonal stationary point: left derivative upper="<<report_upper(dl.upper())
             <<", right derivative lower="<<report_lower(dr.lower())
             <<", value lower="<<report_lower(root_value)<<"\n"
             <<"former equality-point lower="<<report_lower(zstar.v.lower())<<"\n";
    return ok;
}

struct TRBox { double ta,tb,ra,rb; int dep; };
bool local_contains(const TRBox&b) {
    I dmax=I(b.tb)*I(b.rb);
    return b.ta>=LOCAL_T_LO_IN && b.tb<=LOCAL_T_HI_IN
        && dmax.upper()<=LOCAL_D_HI_IN;
}
bool combined_mean_value(const TRBox&b,double&cert) {
    const double tc=0.5*(b.ta+b.tb),rc=0.5*(b.ra+b.rb);
    const double ht=std::nextafter(std::max(tc-b.ta,b.tb-tc),INFINITY);
    const double hr=std::nextafter(std::max(rc-b.ra,b.rb-rc),INFINITY);
    I T(b.ta,b.tb),R(b.ra,b.rb),U=T*(I(1.0)+R),V=T*(I(1.0)-R);
    auto g2=Gall(I(2.0)*T),gu=Gall(U),gv=Gall(V),cu=Call(U),cv=Call(V);
    I one=I(1.0)-THETA;
    I Tc(tc),Rc(rc),Uc=Tc*(I(1.0)+Rc),Vc=Tc*(I(1.0)-Rc);
    I ft0=MTARGET*one*(I(2.0)*ALPHA*g2.p+BETA*((I(1.0)+R)*cu.p*cv.v+(I(1.0)-R)*cu.v*cv.p))
        -HALF*((I(1.0)+R)*gu.p+(I(1.0)-R)*gv.p);
    I fr0=MTARGET*T*one*BETA*(cu.p*cv.v-cu.v*cv.p)-HALF*T*(gu.p-gv.p);
    I fc0=MTARGET*one*(ALPHA*G(I(2.0)*Tc)+BETA*Call(Uc).v*Call(Vc).v)
        -HALF*(G(Uc)+G(Vc));
    I lower0=fc0-iabs(ft0)*I(ht)-iabs(fr0)*I(hr);
    if(lower0.lower()>0){cert=lower0.lower();return true;}

    auto pu=Pall(U),pv=Pall(V); if(!pu.ok||!pv.ok) return false;
    I ft=ft0+MTARGET*THETA*((I(1.0)+R)*pu.p*pv.v+(I(1.0)-R)*pu.v*pv.p);
    I fr=fr0+MTARGET*T*THETA*(pu.p*pv.v-pu.v*pv.p);
    auto puc=Pall(Uc),pvc=Pall(Vc); if(!puc.ok||!pvc.ok) return false;
    I fc=fc0+MTARGET*THETA*puc.v*pvc.v;
    I lower=fc-iabs(ft)*I(ht)-iabs(fr)*I(hr);
    cert=lower.lower(); return cert>0;
}

bool scalar_checks() {
    I c=MTARGET*(I(1.0)-THETA)*ALPHA;
    I t0=decimal_point(0.00001), vb=decimal_point(0.0000001), tailt=decimal_point(8.0);
    I phi=I(2.0)*c*(I(1.0)-ilog(I(2.0)*t0))-I(1.0)
          -iexp(t0)*(-ilog(t0)+t0);
    I boundary=(c-HALF)*G(I(2.0)*t0)-HALF*G(vb);
    I tail=c*(I(2.0)*tailt+I(1.0)-iexp(-I(2.0)*tailt))-(tailt+I(1.0));
    std::cout<<"c_epsilon lower="<<report_lower(c.lower())<<"\n"
             <<"small-t scalar lower="<<report_lower(phi.lower())<<"\n"
             <<"boundary scalar lower="<<report_lower(boundary.lower())<<"\n"
             <<"tail scalar lower="<<report_lower(tail.lower())<<"\n";
    return c.lower()>0.5 && phi.lower()>0 && boundary.lower()>0 && tail.lower()>0;
}

bool certify_combined() {
    // Repeat the scalar check that justifies removing the strip v <= 1e-7.
    // This makes the Proposition A.6 certificate self-contained even though
    // the same numerical inequality is also recorded in Proposition A.4.
    I cc0=MTARGET*(I(1.0)-THETA)*ALPHA;
    I tcorner=decimal_point(0.00001), vcorner=decimal_point(0.0000001);
    I boundary0=(cc0-HALF)*G(I(2.0)*tcorner)-HALF*G(vcorner);
    std::cout<<"A.6 boundary scalar lower="<<report_lower(boundary0.lower())<<"\n";
    if(!(cc0.lower()>0.5 && boundary0.lower()>0)) return false;

    const double T0=down(0.00001),T1=up(8.0),VB=down(0.0000001);
    std::vector<double> tb={T0,0.0001,0.001,0.01,0.05,0.1,0.2,0.3,0.34,0.38,
                            0.45,0.7,1.0,1.5,2.5,4.0,6.0,T1};
    std::vector<TRBox> stack;
    for(size_t i=0;i+1<tb.size();++i) {
        const double a=tb[i],b=tb[i+1];
        const int nr=20;
        I rmaxI=I(1.0)-I(VB)/I(b);
        double rmax=std::min(1.0,rmaxI.upper());
        for(int j=0;j<nr;++j) {
            double ra=(j==0)?0.0:rmax*j/nr;
            double rb=(j+1==nr)?rmax:rmax*(j+1)/nr;
            stack.push_back({a,b,ra,rb,0});
        }
    }
    long long leaves=0,local=0,splits=0; int maxd=0; double mincert=INFINITY;
    TRBox minbox{0,0,0,0,0};
    std::vector<TRBox> bad;
    while(!stack.empty()) {
        TRBox b=stack.back();stack.pop_back();maxd=std::max(maxd,b.dep);
        // Discard boxes wholly covered by the analytic v<=VB estimate.  A box
        // straddling that line is harmless: its interval evaluation proves a
        // slightly larger set and avoids following the curved boundary in
        // (t,xi)-coordinates, where xi is the normalized asymmetry.
        I vhiI=I(b.tb)*(I(1.0)-I(b.ra));
        double vhi=vhiI.upper();
        if(vhi<=VB) continue;
        if(local_contains(b)) {++local;continue;}
        // The same analytic boundary estimate applies with the actual upper
        // endpoint of v.  This removes a wide curved strip near xi=1 without
        // forcing the rectangular subdivision to trace its boundary.
        I cc=MTARGET*(I(1.0)-THETA)*ALPHA;
        I bd=(cc-HALF)*G(I(2.0)*I(b.ta))-HALF*G(I(vhi));
        if(bd.lower()>0) continue;
        double cert=0;
        if(combined_mean_value(b,cert)) {
            if(cert<mincert) {mincert=cert;minbox=b;}
            ++leaves;continue;
        }
        if(b.dep>=38) {bad.push_back(b);if(bad.size()>=20)break;continue;}
        double wt=b.tb-b.ta, wr=(b.rb-b.ra)*b.tb;
        if(wt>=wr) {
            double m=0.5*(b.ta+b.tb);
            stack.push_back({b.ta,m,b.ra,b.rb,b.dep+1});
            stack.push_back({m,b.tb,b.ra,b.rb,b.dep+1});
        } else {
            double m=0.5*(b.ra+b.rb);
            stack.push_back({b.ta,b.tb,b.ra,m,b.dep+1});
            stack.push_back({b.ta,b.tb,m,b.rb,b.dep+1});
        }
        ++splits;
        if(splits%100000==0)
            std::cerr<<"combined progress splits="<<splits<<", stack="<<stack.size()
                     <<", current t=["<<b.ta<<","<<b.tb<<"], xi=["<<b.ra<<","<<b.rb<<"]\n";
    }
    std::cout<<"three-construction inequality: leaves="<<leaves<<", splits="<<splits
             <<", local="<<local<<", maxdepth="<<maxd<<", unresolved="<<bad.size()
             <<", min certified lower="<<report_lower(mincert)<<"\n"
             <<"minimum box: t=["<<minbox.ta<<","<<minbox.tb<<"], xi=["
             <<minbox.ra<<","<<minbox.rb<<"]\n";
    if(!bad.empty()) for(const auto&b:bad)
        std::cout<<"  bad t=["<<b.ta<<","<<b.tb<<"], xi=["<<b.ra<<","<<b.rb<<"]\n";
    return bad.empty();
}

bool verify_constants() {
    auto poly=[](const I&x){return x*x*x*x-I(2.0)*x*x*x+I(3.0)*x*x-I(1.0);};
    I pl=poly(decimal_point(0.69078759392498)),pu=poly(decimal_point(0.69078759392500));
    I dp=I(4.0)*XSTAR*XSTAR*XSTAR-I(6.0)*XSTAR*XSTAR+I(6.0)*XSTAR;
    I A=Hmono(XSTAR*XSTAR), Ap=I(2.0)*XSTAR*Hp(XSTAR*XSTAR);
    I Bp=Hp(gx(XSTAR))*gpx(XSTAR);
    I mcalc=XSTAR*Hmono(XSTAR)/A;
    I betacalc=(A*(I(1.0)/XSTAR+Hp(XSTAR)/Hmono(XSTAR))-Ap)/(Bp-Ap);
    std::cout<<std::setprecision(17)<<"MPFR version="<<mpfr_get_version()
             <<", precision="<<MPFR_PREC<<" bits\n"
             <<"m_* derived=["<<report_lower(mcalc.lower())<<","<<report_upper(mcalc.upper())<<"]\n"
             <<"beta derived=["<<report_lower(betacalc.lower())<<","<<report_upper(betacalc.upper())<<"]\n"
             <<"m_cert interval=["<<report_lower(MTARGET.lower())<<","<<report_upper(MTARGET.upper())<<"]\n";
    std::cout<<"poly(left)=["<<report_lower(pl.lower())<<","<<report_upper(pl.upper())<<"], poly(right)=["
             <<report_lower(pu.lower())<<","<<report_upper(pu.upper())<<"]\n";
    bool c1=pl.upper()<0,c2=pu.lower()>0,c3=dp.lower()>0;
    bool c4=mcalc.lower()>=MSTAR.lower() && mcalc.upper()<=MSTAR.upper();
    bool c5=betacalc.lower()>=BETA.lower() && betacalc.upper()<=BETA.upper();
    std::cout<<"constant checks: root-left="<<c1<<", root-right="<<c2
             <<", derivative="<<c3<<", m-containment="<<c4
             <<", beta-containment="<<c5<<"\n";
    return c1&&c2&&c3&&c4&&c5;
}

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

// MPFR is used for all transcendental endpoint evaluations.  The elementary
// interval type below rounds every algebraic endpoint outward by one binary64
// step after an IEEE-754 operation. Compile with -fno-fast-math.
#if __has_include(<mpfr.h>)
#  include <mpfr.h>
#else
// Fallback ABI declarations for environments that have the MPFR runtime but
// not the development header.  Normal reproducibility builds should install
// libmpfr-dev and use the official <mpfr.h> branch above.
extern "C" {
typedef long int mpfr_prec_t;
typedef long int mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum { MPFR_RNDN=0, MPFR_RNDZ=1, MPFR_RNDU=2, MPFR_RNDD=3, MPFR_RNDA=4, MPFR_RNDF=5, MPFR_RNDNA=-1 } mpfr_rnd_t;
typedef struct { mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp; mp_limb_t *_mpfr_d; } __mpfr_struct;
typedef __mpfr_struct mpfr_t[1]; typedef __mpfr_struct *mpfr_ptr; typedef const __mpfr_struct *mpfr_srcptr;
void mpfr_init2(mpfr_ptr, mpfr_prec_t); void mpfr_clear(mpfr_ptr);
int mpfr_set_d(mpfr_ptr,double,mpfr_rnd_t); int mpfr_set_str(mpfr_ptr,const char*,int,mpfr_rnd_t); int mpfr_set_ui(mpfr_ptr,unsigned long,mpfr_rnd_t);
int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_neg(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr,mpfr_srcptr); double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t); const char* mpfr_get_version(void);
}
#endif

double down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
double up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
struct I {
    double lo,hi;
    I():lo(0),hi(0){}
    I(double x):lo(x),hi(x){}
    I(double a,double b):lo(a),hi(b){if(a>b)throw std::runtime_error("bad interval");}
    double lower()const{return lo;}
    double upper()const{return hi;}
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
static constexpr mpfr_prec_t MPFR_PREC=192;

enum class TF { LOG, EXP, SQRT };
struct MPFRWorkspace {
    mpfr_t a,b;
    MPFRWorkspace(){ mpfr_init2(a,MPFR_PREC); mpfr_init2(b,MPFR_PREC); }
    ~MPFRWorkspace(){ mpfr_clear(a); mpfr_clear(b); }
};
double mpfr_eval(double x, TF f, mpfr_rnd_t rnd) {
    static thread_local MPFRWorkspace w;
    mpfr_set_d(w.a,x,MPFR_RNDN);
    if(f==TF::LOG) mpfr_log(w.b,w.a,rnd); else if(f==TF::EXP) mpfr_exp(w.b,w.a,rnd); else mpfr_sqrt(w.b,w.a,rnd);
    return mpfr_get_d(w.b,rnd);
}
I ilog(const I&x){ if(x.lower()<=0) throw std::runtime_error("log nonpositive"); return I(mpfr_eval(x.lower(),TF::LOG,MPFR_RNDD),mpfr_eval(x.upper(),TF::LOG,MPFR_RNDU)); }
I iexp(const I&x){ return I(mpfr_eval(x.lower(),TF::EXP,MPFR_RNDD),mpfr_eval(x.upper(),TF::EXP,MPFR_RNDU)); }
I isqrt(const I&x){ if(x.lower()<0) throw std::runtime_error("sqrt negative"); return I(mpfr_eval(x.lower(),TF::SQRT,MPFR_RNDD),mpfr_eval(x.upper(),TF::SQRT,MPFR_RNDU)); }

const I MSTAR=decimal_interval(0.61729091208117,0.61729091208136);
const I BETA=decimal_interval(0.1000525598619,0.1000525598639);
const I ALPHA=I(1.0)-BETA;
const I TSTAR_CENTER=decimal_point(0.36992289186771087);
const I LOCAL_RADIUS=decimal_point(0.06);
const I LOCAL_T_LO_ENC=TSTAR_CENTER-LOCAL_RADIUS;
const I LOCAL_T_HI_ENC=TSTAR_CENTER+LOCAL_RADIUS;
const double LOCAL_T_LO_OUT=LOCAL_T_LO_ENC.lower();
const double LOCAL_T_HI_OUT=LOCAL_T_HI_ENC.upper();
const double LOCAL_D_HI_OUT=LOCAL_RADIUS.upper();
const double LOCAL_T_LO_IN=LOCAL_T_LO_ENC.upper();
const double LOCAL_T_HI_IN=LOCAL_T_HI_ENC.lower();
const double LOCAL_D_HI_IN=LOCAL_RADIUS.lower();

I H(const I&z){return -z*ilog(z)-(I(1.0)-z)*ilog(I(1.0)-z);} I Hp(const I&z){return ilog((I(1.0)-z)/z);}
// g is strictly increasing on [0,1].  Evaluate the polynomial on point
// intervals at the two endpoints, with directed algebraic rounding.
I gxraw(const I&x){return x*x*(I(1.0)+(I(1.0)-x)*(I(1.0)-x));}
I gx(const I&x){I a=gxraw(I(x.lower())), b=gxraw(I(x.upper())); return I(a.lower(),b.upper());}
I gpx(const I&x){return I(4.0)*x-I(6.0)*x*x+I(4.0)*x*x*x;} I gppx(const I&x){return I(4.0)-I(12.0)*x+I(12.0)*x*x;}
I G(const I&u){I x=iexp(-u);return H(x)/x;} I C(const I&u){I x=iexp(-u);return isqrt(H(gx(x)))/x;}
I Ftr(const I&T,const I&R){I U=T*(I(1.0)+R),V=T*(I(1.0)-R);return MSTAR*(ALPHA*G(I(2.0)*T)+BETA*C(U)*C(V))-I(0.5)*(G(U)+G(V));}
bool contained(const I&a,const I&b){return a.lower()>=b.lower()&&a.upper()<=b.upper();}

bool verify_constants(){
    const double xl=0.69078759392498,xu=0.69078759392500; auto poly=[](const I&x){return x*x*x*x-I(2.0)*x*x*x+I(3.0)*x*x-I(1.0);};
    I XL=decimal_point(xl),XU=decimal_point(xu),X(XL.lower(),XU.upper());
    I pl=poly(XL),pu=poly(XU); I dp=I(4.0)*X*X*X-I(6.0)*X*X+I(6.0)*X;
    if(!(pl.upper()<0&&pu.lower()>0&&dp.lower()>0))return false;
    I A=H(X*X),Ap=I(2.0)*X*Hp(X*X),Bp=Hp(gx(X))*gpx(X); I mcalc=X*H(X)/A; I betacalc=(A*(I(1.0)/X+Hp(X)/H(X))-Ap)/(Bp-Ap);
    I tcalc=-ilog(X), local_t(LOCAL_T_LO_OUT,LOCAL_T_HI_OUT);
    std::cout<<std::setprecision(17)<<"MPFR version   = "<<mpfr_get_version()<<", transcendental precision="<<MPFR_PREC<<" bits\n";
    std::cout<<"x_* bracket = ["<<xl<<", "<<xu<<"]\nm_* derived   = ["<<report_lower(mcalc.lower())<<", "<<report_upper(mcalc.upper())<<"]\nbeta derived  = ["<<report_lower(betacalc.lower())<<", "<<report_upper(betacalc.upper())<<"]\nt_* derived   = ["<<report_lower(tcalc.lower())<<", "<<report_upper(tcalc.upper())<<"]\n";
    return contained(mcalc,MSTAR)&&contained(betacalc,BETA)&&contained(tcalc,local_t);
}

// Direct high-precision check of the displayed decimal lower bound for c_L.
// The binary64 interval used elsewhere is intentionally wider; this scalar
// check uses MPFR directly so that the last displayed digits are certified.
bool certify_cL_digits(){
    static constexpr mpfr_prec_t P=256;
    mpfr_t x,one,two,three,half,thr,a,b,c,poly,z,Hx,Hx2,num,m,cL,tmp,logz,q,logq,p1,p2;
    mpfr_ptr vars[]={x,one,two,three,half,thr,a,b,c,poly,z,Hx,Hx2,num,m,cL,tmp,logz,q,logq,p1,p2};
    for(auto v:vars) mpfr_init2(v,P);
    auto clear_all=[&](){for(auto v:vars) mpfr_clear(v);};
    mpfr_set_d(x,0.6907875939249879,MPFR_RNDN); // exact binary64 number, below x_*
    mpfr_set_ui(one,1,MPFR_RNDN); mpfr_set_ui(two,2,MPFR_RNDN); mpfr_set_ui(three,3,MPFR_RNDN);
    mpfr_set_d(half,0.5,MPFR_RNDN);
    mpfr_set_str(thr,"0.382709087918735",10,MPFR_RNDU);

    // Upward Horner evaluation of p(x)=((x-2)x+3)x^2-1.
    mpfr_sub(a,x,two,MPFR_RNDU);       // a<0, upper bound
    mpfr_mul(b,a,x,MPFR_RNDU); mpfr_add(b,b,three,MPFR_RNDU); // b>0, upper bound
    mpfr_mul(c,b,x,MPFR_RNDU); mpfr_mul(c,c,x,MPFR_RNDU);
    mpfr_sub(poly,c,one,MPFR_RNDU);
    mpfr_set_ui(tmp,0,MPFR_RNDN);
    bool root_above = mpfr_cmp(poly,tmp)<0;

    // c_L is increasing on the broad x_* bracket used by verify_constants().
    I X(0.69078759392498,0.69078759392500);
    I B=I(1.0)/X+Hp(X)/H(X)-I(2.0)*X*Hp(X*X)/H(X*X);
    bool increasing = B.upper()<0; // m'/m=B, hence c_L'=-m'>0.

    // Upper bound H(x): for a positive point z, evaluating z log z
    // downward and then negating gives an upper bound for -z log z.
    auto entropy_upper=[&](mpfr_ptr out, mpfr_srcptr zz){
        mpfr_log(logz,zz,MPFR_RNDD); mpfr_mul(p1,zz,logz,MPFR_RNDD); mpfr_neg(p1,p1,MPFR_RNDN);
        mpfr_sub(q,one,zz,MPFR_RNDN); // exact by Sterbenz here
        mpfr_log(logq,q,MPFR_RNDD); mpfr_mul(p2,q,logq,MPFR_RNDD); mpfr_neg(p2,p2,MPFR_RNDN);
        mpfr_add(out,p1,p2,MPFR_RNDU);
    };
    // Lower bound H(z) at an exact positive point z.
    auto entropy_lower=[&](mpfr_ptr out, mpfr_srcptr zz){
        mpfr_log(logz,zz,MPFR_RNDU); mpfr_mul(p1,zz,logz,MPFR_RNDU); mpfr_neg(p1,p1,MPFR_RNDN);
        mpfr_sub(q,one,zz,MPFR_RNDN);
        mpfr_log(logq,q,MPFR_RNDU); mpfr_mul(p2,q,logq,MPFR_RNDU); mpfr_neg(p2,p2,MPFR_RNDN);
        mpfr_add(out,p1,p2,MPFR_RNDD);
    };
    entropy_upper(Hx,x);
    // z is a lower bound for x^2. Since z<1/2 and H is increasing on (0,1/2),
    // H(z) is a lower bound for H(x^2).
    mpfr_mul(z,x,x,MPFR_RNDD);
    bool z_below_half=mpfr_cmp(z,half)<0;
    entropy_lower(Hx2,z);
    mpfr_mul(num,x,Hx,MPFR_RNDU);
    mpfr_div(m,num,Hx2,MPFR_RNDU);
    mpfr_sub(cL,one,m,MPFR_RNDD);
    bool digits = mpfr_cmp(cL,thr)>0;
    std::cout<<"direct c_L > 0.382709087918735 : "<<(digits?"yes":"no")<<"\n";
    bool ok=root_above&&increasing&&z_below_half&&digits;
    clear_all(); return ok;
}

struct Triple{I v,p,pp;};
Triple Gall(const I&u){I x=iexp(-u),hh=H(x),hp=Hp(x);I v=hh/x,p=v-hp,pp=p-I(1.0)/(I(1.0)-x);return{v,p,pp};}
Triple Call(const I&u){I x=iexp(-u),z=gx(x),J=H(z),hp=Hp(z),hpp=-I(1.0)/(z*(I(1.0)-z));I z1=-x*gpx(x),z2=x*gpx(x)+x*x*gppx(x);I j1=hp*z1,j2=hpp*z1*z1+hp*z2;I f1=I(1.0)+I(0.5)*j1/J,f2=I(0.5)*(j2/J-(j1/J)*(j1/J));I c=isqrt(J)/x;return{c,c*f1,c*(f2+f1*f1)};}

bool certify_local_hessian(){const int N=40;auto tg=partition(LOCAL_T_LO_OUT,LOCAL_T_HI_OUT,N),dg=partition(0.0,LOCAL_D_HI_OUT,N);double min_tt=INFINITY,min_dd=INFINITY,max_td=0,min_det=INFINITY;bool ok=true;for(int i=0;i<N;i++)for(int j=0;j<N;j++){I T(tg[i],tg[i+1]),D(dg[j],dg[j+1]),U=T+D,V=T-D;auto g2=Gall(I(2.0)*T),gu=Gall(U),gv=Gall(V);auto cu=Call(U),cv=Call(V);I ftt=I(4.0)*MSTAR*ALPHA*g2.pp+MSTAR*BETA*(cu.pp*cv.v+I(2.0)*cu.p*cv.p+cu.v*cv.pp)-I(0.5)*(gu.pp+gv.pp);I fdd=MSTAR*BETA*(cu.pp*cv.v-I(2.0)*cu.p*cv.p+cu.v*cv.pp)-I(0.5)*(gu.pp+gv.pp);I ftd=MSTAR*BETA*(cu.pp*cv.v-cu.v*cv.pp)-I(0.5)*(gu.pp-gv.pp);double lt=ftt.lower(),ld=fdd.lower(),mt=std::max(std::abs(ftd.lower()),std::abs(ftd.upper())); I detI=I(lt)*I(ld)-I(mt)*I(mt); double det=detI.lower();min_tt=std::min(min_tt,lt);min_dd=std::min(min_dd,ld);max_td=std::max(max_td,mt);min_det=std::min(min_det,det);if(lt<=0||ld<=0||det<=0)ok=false;}std::cout<<"local Hessian: min F_tt="<<report_lower(min_tt)<<", min F_dd="<<report_lower(min_dd)<<", max |F_td|="<<report_upper(max_td)<<", min determinant="<<report_lower(min_det)<<"\n";return ok;}
struct Box{double ta,tb,ra,rb;int dep;};
I iabs(const I&x){
    if(x.lower()>=0.0) return x;
    if(x.upper()<=0.0) return -x;
    return I(0.0,std::max(-x.lower(),x.upper()));
}

bool mean_value_positive(const Box& b, double &cert_lower){
    const double tc=0.5*(b.ta+b.tb), rc=0.5*(b.ra+b.rb);
    const double ht=std::nextafter(std::max(tc-b.ta,b.tb-tc),std::numeric_limits<double>::infinity());
    const double hr=std::nextafter(std::max(rc-b.ra,b.rb-rc),std::numeric_limits<double>::infinity());
    I fc=Ftr(I(tc),I(rc));
    I T(b.ta,b.tb), Rr(b.ra,b.rb), U=T*(I(1.0)+Rr), V=T*(I(1.0)-Rr);
    auto g2=Gall(I(2.0)*T), gu=Gall(U), gv=Gall(V);
    auto cu=Call(U), cv=Call(V);
    I ft=I(2.0)*MSTAR*ALPHA*g2.p
        +MSTAR*BETA*((I(1.0)+Rr)*cu.p*cv.v+(I(1.0)-Rr)*cu.v*cv.p)
        -I(0.5)*((I(1.0)+Rr)*gu.p+(I(1.0)-Rr)*gv.p);
    I fr=MSTAR*BETA*T*(cu.p*cv.v-cu.v*cv.p)-I(0.5)*T*(gu.p-gv.p);
    I lower=fc-iabs(ft)*I(ht)-iabs(fr)*I(hr);
    cert_lower=lower.lower();
    return cert_lower>0.0;
}

bool local_contains(const Box&b){I T(b.ta,b.tb),R(b.ra,b.rb),D=T*R;return b.ta>=LOCAL_T_LO_IN&&b.tb<=LOCAL_T_HI_IN&&b.ra>=0.0&&D.upper()<=LOCAL_D_HI_IN;}
bool certify_region(double T0_paper,double T1_paper,double boundary_v_paper,
                    int nbands,int rbins,int maxdepth){
    const double T0=down(T0_paper),T1=up(T1_paper),boundary_v=down(boundary_v_paper);
    auto tg=partition(T0,T1,nbands);
    std::vector<Box> stack;
    stack.reserve(800000);
    for(int i=0;i<nbands;i++){
        double ta=tg[i],tb=tg[i+1];
        I rmaxI=I(1.0)-I(boundary_v)/I(tb);
        double rmax=rmaxI.upper();
        if(rmax<=0) continue;
        rmax=std::min(1.0,rmax);
        auto rg=partition(0.0,rmax,rbins);
        for(int j=0;j<rbins;j++) stack.push_back({ta,tb,rg[j],rg[j+1],0});
    }
    long long leaves=0,local=0,splits=0;
    int maxd=0;
    double min_leaf=INFINITY;
    std::vector<Box> unresolved;
    while(!stack.empty()){
        Box b=stack.back();
        stack.pop_back();
        maxd=std::max(maxd,b.dep);
        if(local_contains(b)){local++;continue;}
        double cert=0.0;
        if(mean_value_positive(b,cert)){
            min_leaf=std::min(min_leaf,cert);
            leaves++;
            continue;
        }
        if(b.dep>=maxdepth){
            unresolved.push_back(b);
            if(unresolved.size()>=20) break;
            continue;
        }
        double wt=b.tb-b.ta,wd=b.tb*b.rb-b.ta*b.ra;
        if(wt>=wd){
            double mid=0.5*(b.ta+b.tb);
            stack.push_back({b.ta,mid,b.ra,b.rb,b.dep+1});
            stack.push_back({mid,b.tb,b.ra,b.rb,b.dep+1});
        }else{
            double mid=0.5*(b.ra+b.rb);
            stack.push_back({b.ta,b.tb,b.ra,mid,b.dep+1});
            stack.push_back({b.ta,b.tb,mid,b.rb,b.dep+1});
        }
        splits++;
    }
    std::cout<<std::setprecision(12)<<"compact ["<<T0_paper<<","<<T1_paper<<"], v>="<<boundary_v_paper
             <<std::setprecision(17)<<": leaves="<<leaves<<", local="<<local<<", splits="<<splits
             <<", maxdepth="<<maxd<<", unresolved="<<unresolved.size()
             <<", min certified lower="<<report_lower(min_leaf)<<"\n";
    return unresolved.empty();
}
bool scalar_checks(){I c=MSTAR*ALPHA,t0=decimal_point(0.002),eps=decimal_point(0.00005),T1=decimal_point(1.4),vb=decimal_point(0.1),Ttail=decimal_point(4.1);I phi=I(2.0)*c*(I(1.0)-ilog(I(2.0)*t0))-I(1.0)-iexp(t0)*(-ilog(t0)+t0);I bd1=(c-I(0.5))*G(I(2.0)*t0)-I(0.5)*G(eps);I bd2=(c-I(0.5))*G(I(2.0)*T1)-I(0.5)*G(vb);I tail=c*(I(2.0)*Ttail+I(1.0)-iexp(-I(2.0)*Ttail))-(Ttail+I(1.0));std::cout<<"small-t eta(0.002) lower = "<<report_lower(phi.lower())<<"\nboundary v<=5e-5 lower    = "<<report_lower(bd1.lower())<<"\nboundary t>=1.4,v<=.1 lower= "<<report_lower(bd2.lower())<<"\ntail t>=4.1 lower          = "<<report_lower(tail.lower())<<"\n";return phi.lower()>0&&bd1.lower()>0&&bd2.lower()>0&&tail.lower()>0;}
int main(){std::cout<<std::setprecision(17);bool ok=true;ok&=verify_constants();ok&=certify_cL_digits();ok&=scalar_checks();ok&=certify_local_hessian();ok&=certify_region(0.002,0.1,0.00005,16,8,22);ok&=certify_region(0.1,0.34,0.00005,20,8,24);ok&=certify_region(0.34,0.40,0.00005,16,8,22);ok&=certify_region(0.40,1.4,0.00005,24,8,24);ok&=certify_region(1.4,4.1,0.1,32,10,22);std::cout<<(ok?"PROPOSITION A.1 CERTIFICATE PASSED\n":"PROPOSITION A.1 CERTIFICATE FAILED\n");return ok?0:1;}

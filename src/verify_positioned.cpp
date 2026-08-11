#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

static constexpr int Q=6, G=4000;
using Range=std::pair<int,int>;
struct Box { std::array<Range,Q> x; uint8_t depth; };
inline double down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
inline long double ldown(long double x){return std::nextafterl(x,-std::numeric_limits<long double>::infinity());}
inline long double lup(long double x){return std::nextafterl(x,std::numeric_limits<long double>::infinity());}

struct I {long double lo,hi;};
inline I point(long double x){return {x,x};}
inline I add(I a,I b){return {ldown(a.lo+b.lo),lup(a.hi+b.hi)};}
inline I neg(I a){return {-a.hi,-a.lo};}
inline I sub(I a,I b){return add(a,neg(b));}
inline I mul(I a,I b){
    long double p[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    return {ldown(*std::min_element(p,p+4)),lup(*std::max_element(p,p+4))};
}
inline I sqr(I a){
    long double hi=std::max(a.lo*a.lo,a.hi*a.hi);
    long double lo=(a.lo<=0&&a.hi>=0)?0:std::min(a.lo*a.lo,a.hi*a.hi);
    return {lo==0?0:ldown(lo),lup(hi)};
}
inline I divpos(I a,I b){ // b.lo > 0
    I inv={ldown(1.0L/b.hi),lup(1.0L/b.lo)};return mul(a,inv);
}
// std::fabs is overloaded for long double; std::fabsl is not surfaced by every
// <cmath> (libstdc++ omits it). Same operation, portable spelling.
inline long double absup(I a){return lup(std::max(std::fabs(a.lo),std::fabs(a.hi)));}

struct RM {
    int n;std::vector<std::vector<double>> lev;
    explicit RM(std::vector<double> v):n((int)v.size()){
        lev.push_back(std::move(v));int w=1;
        while(2*w<=n){const auto&p=lev.back();int half=w;w*=2;std::vector<double> row(n-w+1);for(int i=0;i<=n-w;i++)row[i]=std::min(p[i],p[i+half]);lev.push_back(std::move(row));}
    }
    inline double query(int l,int r)const{if(r>=n)return -std::numeric_limits<double>::infinity();unsigned len=(unsigned)(r-l+1);int k=31-__builtin_clz(len),w=1<<k;const auto&row=lev[k];return std::min(row[l],row[r-w+1]);}
};
static std::vector<double> readbin(const char*path){std::ifstream f(path,std::ios::binary|std::ios::ate);if(!f)throw std::runtime_error(path);size_t bytes=f.tellg();f.seekg(0);std::vector<double> a(bytes/8);f.read((char*)a.data(),bytes);return a;}
struct PairW{int i,j;double lo,hi;};
struct Stats{uint64_t nodes=0,pruned=0,splits=0,tangent=0,convex=0;int depth=0;bool terminal=false;Box bad{};double low=0;};

int main(int argc,char**argv){
    long long tnum=argc>1?std::stoll(argv[1]):509,tden=argc>2?std::stoll(argv[2]):100000;
    int nthreads=argc>3?std::stoi(argv[3]):5;
    std::string table_dir=argc>4?argv[4]:"tables";
    bool use_tangent=!(argc>5 && std::string(argv[5])=="no-tangent");
    double target=(double)tnum/(double)tden,target_up=up(target);
    auto lower=readbin((table_dir+"/w_lower.bin").c_str());RM wrm(std::move(lower));
    auto second=readbin((table_dir+"/w_second_lower.bin").c_str());RM srm(std::move(second));
    auto wlo=readbin((table_dir+"/w_mid_lower.bin").c_str());
    auto whi=readbin((table_dir+"/w_mid_upper.bin").c_str());
    auto dlo=readbin((table_dir+"/w_prime_mid_lower.bin").c_str());
    auto dhi=readbin((table_dir+"/w_prime_mid_upper.bin").c_str());
    constexpr long long WD=1000000000LL,DP=2300000000LL;
    long long ag[12]={219313883,870819909,325978722,WD,WD,2*WD,330070053,129180091,674021278,0,450616064,0};
    long long nums[7][7]{};
    nums[0][1]=nums[5][6]=ag[0];nums[0][2]=nums[4][6]=ag[1];nums[0][3]=nums[3][6]=ag[2];nums[0][4]=nums[2][6]=ag[3];nums[0][5]=nums[1][6]=ag[4];nums[0][6]=ag[5];nums[1][2]=nums[4][5]=ag[6];nums[1][3]=nums[3][5]=ag[7];nums[1][4]=nums[2][5]=ag[8];nums[1][5]=ag[9];nums[2][3]=nums[3][4]=ag[10];nums[2][4]=ag[11];
    std::vector<PairW> pw;for(int i=0;i<=Q;i++)for(int j=i+1;j<=Q;j++)if(nums[i][j])pw.push_back({i,j,down((double)nums[i][j]/WD),up((double)nums[i][j]/WD)});
    long long bn[Q]={831522,1096590,1071888,1071888,1096590,831522};double bp_lo[Q],bp_hi[Q];for(int i=0;i<Q;i++){bp_lo[i]=down((double)bn[i]/DP);bp_hi[i]=up((double)bn[i]/DP);}
    for(int span=1;span<=Q;span++){long long cap=0;for(int i=0;i+span<=Q;i++)cap+=nums[i][i+span];if(cap!=2*WD){std::cerr<<"bad span capacity "<<span<<"\n";return 2;}}
    long long pressure_num=0;for(auto x:bn)pressure_num+=x;if(pressure_num!=6000000){std::cerr<<"bad pressure total\n";return 2;}
    double min_pressure=*std::min_element(bp_lo,bp_lo+Q);
    int required=(int)std::ceil(target_up*G/min_pressure)+1;
    if(wrm.n<required||srm.n<required){std::cerr<<"table too short: have "<<wrm.n<<" need "<<required<<"\n";return 2;}
    if((int)wlo.size()<2*wrm.n+1||whi.size()!=wlo.size()||dlo.size()!=wlo.size()||dhi.size()!=wlo.size()){std::cerr<<"midpoint table mismatch\n";return 2;}

    std::array<std::vector<Range>,Q> comps;
    for(int c=0;c<Q;c++){
        double adjlo=nums[c][c+1]?down((double)nums[c][c+1]/WD):0;bool active=false;int start=0,prev=-2;
        for(int idx=0;idx<wrm.n;idx++){
            double one=down(bp_lo[c]*(double)idx/G);one=down(one+down(adjlo*wrm.lev[0][idx]));
            if(one<target_up){if(!active||idx>prev+1){if(active)comps[c].push_back({start,prev});start=idx;active=true;}prev=idx;}
        }if(active)comps[c].push_back({start,prev});
    }
    std::vector<Box> initial;Box seed{};std::function<void(int)> rec=[&](int c){if(c==Q){seed.depth=0;initial.push_back(seed);return;}for(auto r:comps[c]){seed.x[c]=r;rec(c+1);}};rec(0);
    std::cout<<"target="<<std::setprecision(17)<<target<<" cells="<<wrm.n<<" threads="<<nthreads<<" components=";for(int i=0;i<Q;i++)std::cout<<comps[i].size()<<(i<5?",":"");std::cout<<" initial="<<initial.size()<<" tangent="<<(use_tangent?"on":"off")<<" required_cells="<<required<<"\n"<<std::flush;

    auto signed_scalar=[&](const PairW&p,double s){double factor=s>=0?p.lo:p.hi;return down(factor*s);};
    auto prove_pd=[&](const Box&box,int*lp,int*hp)->bool{
        double A[Q][Q]{};I AI[Q][Q];for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)AI[i][j]=point(0);
        for(const auto&p:pw){int span=p.j-p.i,l=lp[p.j]-lp[p.i],r=hp[p.j]-hp[p.i]+span-1;double s=srm.query(l,r);if(!std::isfinite(s))return false;double c=signed_scalar(p,s);I ci=point((long double)c);for(int u=p.i;u<p.j;u++)for(int v=p.i;v<p.j;v++){A[u][v]+=c;AI[u][v]=add(AI[u][v],ci);}}
        // Cheap numerical gate.
        double L[Q][Q]{},D[Q]{};
        for(int col=0;col<Q;col++){double piv=A[col][col];for(int k=0;k<col;k++)piv-=L[col][k]*L[col][k]*D[k];if(!(piv>1e-12))return false;D[col]=piv;L[col][col]=1;for(int row=col+1;row<Q;row++){double v=A[row][col];for(int k=0;k<col;k++)v-=L[row][k]*L[col][k]*D[k];L[row][col]=v/piv;}}
        // Rigorous outward-rounded interval LDL.
        I LI[Q][Q],DI[Q];for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)LI[i][j]=point(0);
        for(int col=0;col<Q;col++){
            I piv=AI[col][col];for(int k=0;k<col;k++)piv=sub(piv,mul(sqr(LI[col][k]),DI[k]));if(!(piv.lo>0))return false;DI[col]=piv;LI[col][col]=point(1);
            for(int row=col+1;row<Q;row++){I v=AI[row][col];for(int k=0;k<col;k++)v=sub(v,mul(mul(LI[row][k],LI[col][k]),DI[k]));LI[row][col]=divpos(v,piv);}
        }return true;
    };
    auto tangent_lower=[&](const Box&box,int*lp,int*hp)->long double{
        int midnum[Q];I value=point(0),grad[Q];for(int c=0;c<Q;c++){midnum[c]=box.x[c].first+box.x[c].second+1;I p={(long double)bp_lo[c],(long double)bp_hi[c]};I x={ldown((long double)midnum[c]/(2*G)),lup((long double)midnum[c]/(2*G))};value=add(value,mul(p,x));grad[c]=p;}
        int mp[Q+1]={0};for(int c=0;c<Q;c++)mp[c+1]=mp[c]+midnum[c];
        for(const auto&p:pw){int ix=mp[p.j]-mp[p.i];if(ix<0||ix>=(int)wlo.size())return -std::numeric_limits<long double>::infinity();I a={(long double)p.lo,(long double)p.hi};I w={(long double)wlo[ix],(long double)whi[ix]};I der={(long double)dlo[ix],(long double)dhi[ix]};value=add(value,mul(a,w));I term=mul(a,der);for(int c=p.i;c<p.j;c++)grad[c]=add(grad[c],term);}
        long double out=value.lo;
        for(int c=0;c<Q;c++){long double rad=lup((long double)(box.x[c].second-box.x[c].first+1)/(2*G));out=ldown(out-lup(absup(grad[c])*rad));}
        return out;
    };

    std::atomic<size_t> next{0};std::atomic<bool>stop{false};std::vector<Stats> all(nthreads);auto t0=std::chrono::steady_clock::now();
    auto worker=[&](int tid){Stats s;std::vector<Box> stack;stack.reserve(256);while(!stop.load(std::memory_order_relaxed)){size_t ix=next.fetch_add(1);if(ix>=initial.size())break;stack.clear();stack.push_back(initial[ix]);while(!stack.empty()&&!stop.load(std::memory_order_relaxed)){Box box=stack.back();stack.pop_back();s.nodes++;s.depth=std::max(s.depth,(int)box.depth);int lp[Q+1]={0},hp[Q+1]={0};double pres=0;for(int c=0;c<Q;c++){int lo=box.x[c].first,hi=box.x[c].second;lp[c+1]=lp[c]+lo;hp[c+1]=hp[c]+hi;pres=down(pres+down(bp_lo[c]*(double)lo/G));}if(pres>=target_up){s.pruned++;continue;}double low=pres;for(const auto&p:pw){int span=p.j-p.i,l=lp[p.j]-lp[p.i],r=hp[p.j]-hp[p.i]+span-1;double km=wrm.query(l,r);if(!std::isfinite(km))km=0;low=down(low+down(p.lo*km));}if(low>=target_up){s.pruned++;continue;}
                if(use_tangent && prove_pd(box,lp,hp)){s.convex++;long double tl=tangent_lower(box,lp,hp);if(tl>=(long double)target_up){s.pruned++;s.tangent++;continue;}}
                int k=0,wmax=-1;for(int c=0;c<Q;c++){int w=box.x[c].second-box.x[c].first;if(w>wmax){wmax=w;k=c;}}if(wmax==0){s.terminal=true;s.bad=box;s.low=low;stop.store(true);break;}int l=box.x[k].first,r=box.x[k].second,m=(l+r)/2;Box lo=box,hi=box;lo.depth=hi.depth=box.depth+1;lo.x[k]={l,m};hi.x[k]={m+1,r};stack.push_back(lo);stack.push_back(hi);s.splits++;}}
        all[tid]=s;};
    std::vector<std::thread>th;for(int i=0;i<nthreads;i++)th.emplace_back(worker,i);for(auto&x:th)x.join();double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();Stats z;for(auto&s:all){z.nodes+=s.nodes;z.pruned+=s.pruned;z.splits+=s.splits;z.tangent+=s.tangent;z.convex+=s.convex;z.depth=std::max(z.depth,s.depth);if(s.terminal){z.terminal=true;z.bad=s.bad;z.low=s.low;}}
    if(z.terminal){std::cout<<"VERIFIED=false terminal_low="<<std::setprecision(17)<<z.low<<" box=";for(auto&r:z.bad.x)std::cout<<"["<<r.first<<","<<r.second<<"]";std::cout<<" nodes="<<z.nodes<<"\n";return 1;}
    std::cout<<"VERIFIED=true nodes="<<z.nodes<<" pruned="<<z.pruned<<" splits="<<z.splits<<" convex="<<z.convex<<" tangent="<<z.tangent<<" max_depth="<<z.depth<<" seconds="<<sec<<" rate="<<z.nodes/sec<<"\n";return 0;
}

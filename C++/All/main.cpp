#include <bits/stdc++.h>

#define ll long long
#define db double long
#define pii pair<int,int>
#define pll pair<ll,ll>
#define name "main"
#define inout freopen(name".inp","r", stdin); freopen(name".out","w",stdout);

using namespace std;

const int N=1e4+10;
const int M=log2(N)+10;
const ll MOD=111539786;
const int base=257;
const int INF=1e9;
const ll INFI=1e18;
const db pi=acos(-1.00);
const db eps=1e-6;

int dx4[]={-1, 1, 1, 0, 0};
int dy4[]={0, 0, -1, 1};
int dx8[]={-1, -1, 0, 1, 1, 1, 0, -1};
int dy8[]={0, 1, 1, 1, 0, -1, -1, -1};

template <class T> ll sqr(T a) { return 1LL*a*a; }
template<typename T, typename U> inline void amin(T &x, U y) { if(y < x) x = y; }
template<typename T, typename U> inline void amax(T &x, U y) { if(x < y) x = y; }
template<typename T> inline void read(T &x){
    x = 0; char c;
    while(!isdigit(c = getchar()));
    do
        x = x * 10 + c - '0';
    while(isdigit(c = getchar()));
}
template<typename T> inline void write(T x){
    if (x > 9) write(x / 10);
    putchar(x % 10 + '0');
}

#define Dương Minh Lợi 21120017 - 21CNTN HCMUS

///* List các thuật ghi nhớ *///

/// ___ Toán  ___
/// n^k (Tính n^k)
/// C(k,n) (Tính tổ hợp chập k của n)
/// Sieve (Sàng nguyên tố)
/// Factorize (Phân tích số thành số nguyên tố)
/// Divisors (Đếm số ước)
/// Convex hull (Bao lồi)


/// ___ Quy hoạch động ___
/// Matrix (Nhân ma trận)



/// ___ Đồ thị ___
/// LCA
/// Bridge Cut (Tìm  khớp  cầu)
/// Topo sort


/// ___ Cấu trúc dữ liệu ___
/// Sparse Table
/// Deque
/// Trie

#define Dương Minh Lợi 12 Toán 1 Trường THPT Chuyên Lê Quý Đôn

/// n^k trong log(k)

ll Pow(ll n, ll k){
    if (k==1){
        return n;
    }
    ll res=Pow(n,k/2);
    res*=res;
    res%=MOD;
    if (k&1){
        n%=MOD;
        res*=n;
        res%=MOD;
    }
    return res;
}

/// Giai thừa
/// fac[i] là giai thừa đến i
/// ifac[i] là nghịch đảo modulo
/// ifac[i]=pow(fac[i],mod-2) vì ifac[i]=pow(fac[i],-1) mà mod-2==-1 theo phi hàm euler

ll fac[N+10], ifac[N+10];

void preC(){
    fac[0]=1;
    for (int i=1; i<=N; i++){
        fac[i]=fac[i-1]*i;
        fac[i]%=MOD;
    }
    ifac[N]=Pow(fac[N],MOD-2);
    for (int i=N-1; i>=0; i--){
        ifac[i]=ifac[i+1]*(i+1);
        ifac[i]%=MOD;
    }
}

ll C(ll k, ll n){
    ll res=fac[n];
    res*=ifac[k];
    res%=MOD;
    res*=ifac[n-k];
    res%=MOD;
    return res;
}

/// Sàng nguyên tố
/// prime[i] là số nguyên tố nhỏ nhất mà i chia hết

int prime[N];

void Sieve(){
    for (int i=2; i<=N; i++){
        prime[i]=INF;
    }
    for (int i=2; i*i<=N; i++){
        if (prime[i]==INF){
            prime[i]=i;
            for (int j=i*i; j<=N; j+=i){
                prime[j]=min(prime[j],i);
            }
        }
    }
    for (int i=2; i<=N; i++){
        if (prime[i]==INF){
            prime[i]=i;
        }
    }
}

/// Phân tích thừa số nguyên tố

vector <int> factorize(int n){
    vector <int> res;
    while (n!=1){
        res.push_back(prime[n]);
        n/=prime[n];
    }
    return res;
}

/// Tìm số ước

int divisors(int n){
    vector <int> f=factorize(n);
    int res=1;
    int cnt=1;
    f.push_back(0);
    for (int i=0; i<=f.size()-2; i++){
        cnt++;
        if (f[i]!=f[i+1]){
            res*=cnt;
            cnt=1;
        }
    }
    return res;
}

/// Convex hull

struct point{
    ll x,y;
};

bool cmp(point A, point B){
    return (A.x==B.x) ? (A.y<B.y) : (A.x<B.x);
}

ll cross(point A, point B, point C){
    return A.x*(B.y-C.y)+B.x*(C.y-A.y)+C.x*(A.y-B.y);
}

point con[N], p[N];

void RunConvex_hull(){
    while (true){
        int n;
        cin >> n;
        if (!n){
            break;
        }
        for (int i=1; i<=n; i++){
            cin >> p[i].x >> p[i].y;
        }
        sort(p+1,p+n+1,cmp);
        int k=0;
        for (int i=1; i<=n; i++){
            while (k>=2 && cross(con[k-1],con[k],p[i])<=0){
                k--;
            }
            con[++k]=p[i];
        }
        int t=k+1;
        for (int i=n; i>=1; i--){
            while (k>=t && cross(con[k-1],con[k],p[i])<=0){
                k--;
            }
            con[++k]=p[i];
        }
        k--;
        vector <point> res;
        for (int i=1; i<=k; i++){
            if (!res.size()){
                res.push_back(con[i]);
            } else {
                point last=res.back();
                if (last.x!=con[i].x || last.y!=con[i].y){
                    res.push_back(con[i]);
                }
            }
        }
        cout << res.size() << "\n";
        for (auto v: res){
            cout << v.x << " " << v.y << "\n";
        }
    }
}

/// Matrix
/// Nhân ma trận chỉ áp dụng với phép "+"
/// Một số bài ko pải "+" thì chuyển về "+"
/// Ví dụ như "*" dp[i]=dp[i-1]*dp[i-2] thì nhân ma trận mũ của nó
/// Xây ma trận = số giá trị trước nó

struct Matrix{
    ll x[2][2];
    Matrix(){
        memset(x,0,sizeof(x));
    }
    friend Matrix operator * (const Matrix &A, const Matrix &B){
        Matrix res;
        for (int i=0; i<=1; i++){
            for (int j=0; j<=1; j++){
                for (int k=0; k<=1; k++){
                    res.x[i][j]+=A.x[i][k]*B.x[k][j];
                    res.x[i][j]%=MOD;
                }
            }
        }
        return res;
    }
};

Matrix Pow(Matrix A, ll k){
    if (k==1){
        return A;
    }
    Matrix res=Pow(A,k/2);
    res=res*res;
    if (k&1){
        res=res*A;
    }
    return res;
}

/// Ma trận cơ sở
/// Lấy cnt giá trị đầu của hàm (cnt==số hàng trong ma trận chuyển)
/// Lúc xây ma trận cơ sở xây từ trên xuống dưới chứ ko pải xây từ trái quá phải
/// Ví dụ mày muốn xây dp[3]=10; dp[2]=4; dp[1]=1;
/// 10 10 10
/// 4 4 4
/// 1 1 1

Matrix cs;

/// Bảng ma trận Chuyển
/// 1*dp[i]  =1*dp[i-1]+1*dp[i-2]   1 1 1
/// 1*dp[i-1]=1*dp[i-1]+0*dp[i-2] = 1 1 0
/// 1*dp[i-2]=0*dp[i-1]+1*dp[i-2]   1 0 1

Matrix c;

void SetMatrix(){
    cs.x[0][0]=2; cs.x[0][1]=2;
    cs.x[1][0]=1; cs.x[1][1]=1;

    c.x[0][0]=1; c.x[0][1]=1;
    c.x[1][0]=1; c.x[1][1]=0;
}

/// Xét mấy trường hợp đầu riêng ra vì ko cái mảng chuyển nó chảy âm đó

void runMatrix(){
    int t;
    cin >> t;
    while (t--){
        int n;
        cin >> n;
        if (n==1){
            cout << 1 << "\n";
            continue;
        }
        if (n==2){
            cout << 2 << "\n";
            continue;
        }
        SetMatrix();
        Matrix res=Pow(c,n-2);
        res=res*cs;
        cout << res.x[0][0] << "\n";
    }
}

/// LCA
vector <int> G[N*100];
int depth[N], anc[N][M],n;

void dfs1(int p, int u){
    depth[u]=depth[p]+1;
    anc[u][0]=p;
    for (auto v: G[u]){
        if (p!=v){
            dfs1(u,v);
        }
    }
}

void prelca(){
    for (int j=1; (1 << j) <=n; j++){
        for (int i=1; i<=n; i++){
            anc[i][j]=anc[anc[i][j-1]][j-1];
        }
    }
}

int lca(int u, int v){
    if (depth[u]<depth[v]){
        swap(u,v);
    }
    for (int i=log2(n); i>=0; i--){
        if (depth[u]-(1<<i)>=depth[v]){
            u=anc[u][i];
        }
    }
    if (u==v){
        return u;
    }
    for (int i=log2(n); i>=0; i--){
        if (anc[u][i]!=anc[v][i]){
            u=anc[u][i];
            v=anc[v][i];
        }
    }
    return anc[u][0];
}

void runLCA(){
    cin >> n;
    for (int i=1; i<n; i++){
        int u,v;
        cin >> u >> v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    dfs1(0,1);
    prelca();
    int q;
    cin >> q;
    while (q--){
        int u,v;
        cin >> u >> v;
        cout << lca(u,v) << "\n";
    }
}

/// Bridge, Cut
/// Cái cầu không cần lưu lại 2 đỉnh. Chỉ cần lưu đỉnh sau rùi lấy parent của nó là dc
/// Cầu là u u. Khớp là u v. Cả 2 cái đều là low>=num

int cnt=0, parent[N], num[N], low[N];
int cut[N], bridge[N];

void dfs2(int p, int u){
    parent[u]=p;
    num[u]=low[u]=++cnt;
    int child=0;
    for (auto v: G[u]){
        if (!num[v]){
            child++;
            dfs2(u,v);
            low[u]=min(low[u],low[v]);
            if (low[v]>=num[v]){
                bridge[v]=1;
            }
            if ((low[v]>=num[u] && p) || (child>=2 && !p)){
                cut[u]=1;
            }
        } else if (v!=p){
            low[u]=min(low[u],num[v]);
        }
    }
}

void runBC(){
    int m;
    cin >> n >> m;
    for (int i=1; i<=m; i++){
        int u,v;
        cin >> u >> v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    for (int i=1; i<=n; i++){
        if (!num[i]){
            dfs2(0,i);
        }
    }
    int Bridge=0;
    int Cut=0;
    for (int i=1; i<=n; i++){
        if (bridge[i]){
            Bridge++;
        }
        if (cut[i]){
            Cut++;
            cout << i << " ";
        }
    }
    //cout << Cut << " " << Bridge << "\n";
}

/// Topo sort leexical
/// Đồ thị có hướng

vector <int> Topo;

int deg[N];

void runTopo(){
    int n,m;
    cin >> n >> m;
    for (int i=1; i<=m; i++){
        int u,v;
        cin >> u >> v;
        G[u].push_back(v);
        deg[v]++;
    }
    int root=0;
    priority_queue <int, vector<int>, greater<int>> pq;
    for (int i=1; i<=n; i++){
        if (!deg[i]){
            pq.push(i);
        }
    }
    while (pq.size()){
        int u=pq.top();
        Topo.push_back(u);
        pq.pop();
        for (auto v: G[u]){
            deg[v]--;
            if (!deg[v]){
                pq.push(v);
            }
        }
    }
    if (Topo.size()!=n){
        cout << "Sandro fails.";
    } else {
        for (auto v: Topo){
            cout << v << " ";
        }
    }
}
/// Sparse table
/// chạy query thì  nhớ chạy cái sau là lấy r-(1<<j)+1

int Minv[N][M], a[N];

void preSparse(){
    for (int i=1; i<=n; i++){
        Minv[i][0]=a[i];
    }
    for (int j=1; j<=log2(n); j++){
        for (int i=1; i+(1<<j)-1<=n; i++){
            Minv[i][j]=min(Minv[i][j-1],Minv[i+(1<<(j-1))][j-1]);
        }
    }
}

int Query(int l, int r){
    int len=log2(r-l+1);
    //cout << "? " << r-(1<<j)+1 << "\n";
    return min(Minv[l][len],Minv[r-(1<<len)+1][len]);
}

void runSparse(){
    cin >> n;
    for (int i=1; i<=n; i++){
        cin >> a[i];
    }
    preSparse();
    int q;
    cin >> q;
    while (q--){
        int l,r;
        cin >> l >> r;
        cout << Query(l,r) << "\n";
    }
}

/// Deque
/// Mẹo nhớ ...... Lấy min tịnh tiến thì đằng trc pải bé hơn...
/// Xét tất cả đoạn con có độ dài k, tìm min-max

void runDeque(){
    int n,pos[N];
    cin >> n;
    for (int i=1; i<=n; i++){
        cin >> a[i];
    }
    deque <int> d;
    for (int i=1; i<=n; i++){
        while (d.size() && a[d.back()]>=a[i]){
            d.pop_back();
        }
        /*while (d.size() && d.front()<i){
            d.pop_front();
        }*/
        if (!d.size()){
            pos[i]=0;
        } else {
            pos[i]=d.back();
        }
        d.push_back(i);
    }
    for (int i=1; i<=n; i++){
        cout << pos[i] << " ";
    }
}

/// Trie
/// Lưu ý số Node là N*30

struct Trie{
    int child[30];
    bool isend;
};

Trie node[N*30];
int trie=0;

void update(string s){
    int cur=0;
    for (char c: s){
        int x=(int)c-'a'+1;
        if (!node[cur].child[x]){
            node[cur].child[x]=++trie;
        }
        cur=node[cur].child[x];
    }
    node[cur].isend=1;
}

bool check(string s){
    int cur=0;
    for (char c: s){
        int x=(int)c-'a'+1;
        if (!node[cur].child[x]){
            return 0;
        }
        cur=node[cur].child[x];
    }
    return node[cur].isend;
}

void runTrie(){
    int n;
    cin >> n;
    for (int i=1; i<=n; i++){
        string s;
        cin >> s;
        update(s);
    }
    int q;
    cin >> q;
    while (q--){
        string s;
        cin >> s;
        cout << check(s) << "\n";
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    //inout;
    //runMatrix();
    //runLCA();
    //runSparse();
    //runBC();
    //RunConvex_hull();
    //runDeque();
    //runTrie();
    //runTopo();
    return 0;
}


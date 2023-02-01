#include<bits/stdc++.h>


using namespace std;

typedef pair<int,int> fraction;

int github_try;
int order, fn,n;
vector<fraction> F;
vector<vector<int>> AFT;

fraction compare( int e, int f ){
    int g = __gcd(e,f);
    return {e/g,f/g};
}

void print_farey(){

    for( auto i : F )    cout << i.first << "/" << i.second << ", ";
    cout << endl;
}
void print_AFT(){
    int p = AFT.size(), q = AFT[0].size();
    for( int i = 0; i < p; ++i ){

        for( int j = 0; j < q; ++j )
            cout << AFT[i][j] << " ";
        cout << endl;
    }
}

int Fn( int x, int y, int val = 0 ){
    int X = abs(x-order);
    int Y = y+order;
    if( val )
        return AFT[X][Y] = val;
    return AFT[X][Y];
}
void extend_AFT(){

    for( int a = 0; a <= n; ++a )
        for( int b = -n; b <= -1; ++b )
            Fn(a,b,4*fn-2-Fn(a,-b));

    for( int a = -n; a <= -1; ++a )
        for( int b = -n; b <= -1; ++b )
            Fn(a,b,8*fn-6-Fn(-a,b));

    for( int a = -n; a <= -1; ++a )
        for( int b = 0; b <= n; ++b )
            Fn(a,b,8*fn-6-Fn(-a,b));
}
void gen_AFT(){

    //DEBUG //cout << 28 << endl;
    n = order;
    int a = 0, b = 1, c = 1, d = n, e, f,k, g, rank = 2;
    for( int i = 1; i <= n; ++i ){
        Fn(0,i,1);
    }

    //DEBUG //cout << 35 << endl;

    F.push_back({0,1}), F.push_back({1,n}), Fn(1,n,2);

    //DEBUG //cout << 39 << endl;

    do{
        rank = rank + 1;
        k = (n+b)/d;
        e = k*c-a, f = k*d-b;
        g = __gcd(e,f);
        e = e/g,    f = f/g;

        //DEBUG //cout << rank << " --> " << a<<"/"<<b<<", " << c<<"/"<<d<<", " << e<<"/"<<f << endl;

        F.push_back({e,f}); Fn(e,f,rank);
        for( int i = 1; i <= n/f; ++i )
            Fn(e*i,f*i,rank);
        a = c,  b = d,  c = e,  d = f;
    }
    while( e != 1 || f != 2 );

    //DEBUG //cout << 54 << endl;

    fn = 2*rank-1;

    for( a = 1; a <= n; ++a ){

        for( b = a; b <= 2*a-1; ++b ){
            Fn(a,b,fn+1-Fn(b-a,b));
        }
    }

    //DEBUG //cout << 65 << endl;

    for( a = 1; a <= n; ++a ){

        for( b = 0; b <= a-1; ++b ){
            Fn(a,b,2*fn-Fn(b,a));
        }
    }

    extend_AFT();

}

double Fraction_Difference( double a1, double b1, double a2, double b2 ){
    return (a1/b1)-(a2/b2);
}
bool check( double x, double y, double a, double b, double X, double Y ){
    return x/y <= a/b && a/b <= X/Y;
}
bool LessThan( double x, double y, double X, double Y ){
    return x/y < X/Y;
}
int find_closest_rank( fraction a_b ){
    int a = a_b.first, b = a_b.second;
    int Floor = (a*n)/b, Ceil = (a*n+n-1)/b;
    int g1 = __gcd( Floor, n ), g2 = __gcd( Ceil, n );
    int a1 = Floor/g1, b1 = n/g1;
    int a2 = Ceil/g2, b2 = n/g2;
    double r1 = Fn(a1,b1), r2 = Fn(a2,b2);

    if( r1 == r2 )
        return Fn(a1,b1);

    if( r2-r1 == 1 ){
        if( Fraction_Difference(a,b,a1,b1) < Fraction_Difference(a2,b2,a,b) )
            return Fn(a1,b1);
        else
            return Fn(a2,b2);
    }

    int a1_ = a1, b1_ = b1, a2_ = a2, b2_ = b2;
    do{

        if( LessThan(a,b,a2_,b2_) )
            a2 = a2_, b2 = b2_;
        else
            a1 = a2_, b1 = b2_;

        double y1 = Fraction_Difference(a1,b1,a,b), y2 = Fraction_Difference(a2,b2,a,b);
        double r = ( r1*y2 - r2*y1 )/(y2-y1);
        a1_ = F[floor(r)].first, b1_ = F[floor(r)].second;
        if( LessThan(a1_,b1_,a,b) )
            a2_ = F[floor(r+(double)1)].first, b2_ = F[floor(r+(double)1)].second;
        else
            a2_ = F[floor(r-(double)1)].first, b2_ = F[floor(r-(double)1)].second;
    }
    while( check(a1_,b1_,a,b,a2_,b2_) || check(a2_,b2_,a,b,a1_,b1_) );

    if( Fraction_Difference(a,b,a1_,b1_) < Fraction_Difference(a2_,b2_,a,b) )
        return Fn(a1_,b1_);

    return Fn(a2_,b2_);
}
int main(){

    F.push_back({-1,-1});

    cout << "Order of AFT:";
    cin >> order;

    AFT = vector<vector<int>> (2*order+1,vector<int>(2*order+1));


    gen_AFT();

    cout << endl;

    print_farey();
    print_AFT();

    cout << find_closest_rank({9,16});
    return 0;
}

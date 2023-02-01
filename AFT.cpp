#include<bits/stdc++.h>
#include <numeric>

#define double long double
using namespace std;

typedef pair<int,int> fraction;

int order, fn,n;
vector<fraction> F;
vector<vector<int>> AFT;

int LCM( int a, int b ){
    return (a*b)/__gcd(a,b);
}
int LCM( int a, int b, int c ){
    int l = LCM(a,b);
    return (l*c)/__gcd(l,c);
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

    assert( abs(x) <= n && abs(y) <= n );

    int X = abs(x-order);
    int Y = y+order;
    //DEBUG //cout << X << " " << Y << endl;
    if( val )
        return AFT[X][Y] = val;
    return AFT[X][Y];
}
void extend_AFT(){

    //DEBUG //cout << "extended_AFT()" << endl;
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
        ////DEBUG //cout << "          " << 81 << endl;
        F.push_back({e,f}); Fn(e,f,rank);
        for( int i = 1; i <= n/f; ++i ){
            //DEBUG //cout << 84 << "  ";
            Fn(e*i,f*i,rank);
        }
        a = c,  b = d,  c = e,  d = f;
    }
    while( e != 1 || f != 2 );

    //DEBUG //cout << 54 << endl;

    fn = 2*rank-1;

    for( a = 1; a <= n; ++a ){
        //DEBUG //cout << 96 << "  ";
        for( b = a; (b+1)/2 <= a && a<=b && b <= n; ++b ){
             //DEBUG//cout << a <<"," << b << "  ";
             Fn(a,b,fn+1-Fn(b-a,b));
        }
    }

    //DEBUG //cout << 65 << endl;

    for( a = 1; a <= n; ++a ){

        for( b = 0; b <= a-1; ++b ){
            //DEBUG //cout << 107 << "  ";
            Fn(a,b,2*fn-Fn(b,a));
        }
    }

    extend_AFT();

}

double Fraction_Difference( double a1, double b1, double a2, double b2 ){
    double l = LCM(b1,b2);
    a1 = a1*(l/b1);
    a2 = a2*(l/b2);
    //return (a1/b1)-(a2/b2);
    return (a1-a2)/l;
}
//Performing LCM and making denominator same to avoid precision error
bool check( int x, int y, int a, int b, int X, int Y ){
    int l = LCM(y,b,Y);
    x = x*(l/y);
    a = a*(l/b);
    X = X*(l/X);
    //denominator of all fractions is now same, hence only the numerator need to be compared
    return x <= a && a <= X;
}
bool LessThan( int x, int y, int X, int Y ){
    int l = LCM(y,Y);
    x = x*(l/y);
    X = X*(l/Y);
    //return x/y < X/Y;
    //denominator of all fractions is now same, hence only the numerator need to be compared
    return x < X;
}
bool LessThan( fraction f1, fraction f2 ){
    int x = f1.first, y = f1.second, X = f2.first, Y = f2.second;
    return LessThan(x,y,X,Y);
}
//DEBUG //void test( int a, int b, int a1, int b1, int a1_, int b1_, int a2, int b2, int a2_, int b2_ ){
    //DEBUG //cout << "a/b - " << a << "/" << b << endl;
    //DEBUG //cout << "a1/b1 - " << a1 << "/" << b1 << "   " << "a1'/b1' - " << a1_ << "/" << b1_ << endl;
    //DEBUG //cout << "a2/b2 - " << a2 << "/" << b2 << "   " << "a2'/b2' - " << a2_ << "/" << b2_ << endl;
//DEBUG //}
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

    //DEBUG //cout << "a=" << a << " b=" << b << endl;
    //DEBUG //cout << "Floor=" << Floor << " Ceil=" << Ceil << endl;
    //DEBUG //cout << "g1=" << g1 << " g2=" << g2 << endl;
    //DEBUG //cout << "a1=" << a1 << " b1=" << b1 << endl;
    //DEBUG //cout << "a2=" << a2 << " b2=" << b2 << endl;
    //DEBUG //cout << "r1=" << r1 << " r2=" << r2 << endl;
    int a1_ = a1, b1_ = b1, a2_ = a2, b2_ = b2;
    //DEBUG //cout << "a1_=" << a1_ << " b1_=" << b1_ << endl;
    //DEBUG //cout << "a2_=" << a2_ << " b2_=" << b2_ << endl;
    //DEBUG //cout << endl << endl;
    int cnt = 0;
    do{

        if( LessThan(a,b,a2_,b2_) )
            a2 = a2_, b2 = b2_;
        else
            a1 = a2_, b1 = b2_;

        //DEBUG //cout << "\t\t17" << endl;
        //DEBUG //test( a, b, a1, b1, a1_, b1_, a2, b2, a2_, b2_ );

        double y1 = Fraction_Difference(a1,b1,a,b), y2 = Fraction_Difference(a2,b2,a,b);
        //DEBUG //cout << "\ty1=" << y1 << " y2=" << y2 << endl;
        double r = ( r1*y2 - r2*y1 )/(y2-y1);
        //DEBUG //cout << "r=" << r << endl;
        a1_ = F[floor(r)].first, b1_ = F[floor(r)].second;

        //DEBUG //cout << "\t\t20" << endl;
        //DEBUG //test( a, b, a1, b1, a1_, b1_, a2, b2, a2_, b2_ );

        if( LessThan(a1_,b1_,a,b) ){
            int idx = (int)floor(r+(double)1);//min( (int)F.size()-1,(int)floor(r+(double)1) );
            a2_ = F[idx].first, b2_ = F[idx].second;
        }
        else{
            int idx = (int)floor(r-(double)1);//max(0,(int)floor(r-(double)1));
            a2_ = F[idx].first, b2_ = F[idx].second;
        }

        //DEBUG //cout << "\t\t24" << endl;
        //DEBUG //test( a, b, a1, b1, a1_, b1_, a2, b2, a2_, b2_ );

        int rank1 = Fn(a1_,b1_), rank2 = Fn(a2_,b2_);
        if( abs(rank1-rank2) == 1 ) break;
        //++cnt;
    }
    while( ( check(a1_,b1_,a,b,a2_,b2_) || check(a2_,b2_,a,b,a1_,b1_) ) && cnt < 10 );

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

    //cout << find_closest_rank({78,261});
    cout << find_closest_rank({3,5});
    return 0;
}

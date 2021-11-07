#include <bits/stdc++.h>
using namespace std ;
#define Day (86164.09054)   /// ������������ �����
#define Me (5.9726e24)      /// ����� �����
#define Re (6378100.0)      /// �������������� ������
#define G (6.6743015e-11)   /// �������������� ����������
#define DT (0.001)          /// ��� �� �������
#define TICKER (10000)    /// ��������� ������������� ������� ���������� � ������� (� �����)

#define Omegae (2*M_PI/Day) /// ������� �������� �����
#define Ro (42165188.26463) /// ������ ��������������� ������
#define Vo (Omegae*Ro)      /// ����������� �������� �� ���������������

#define K (0.05)            /// �������� ����� �����

#define ALPHA (59.0)        /// ���� �������� ������� �������� ������� �����
#define T3    (35.335)      /// ������������ ������� ��������� ��������� 3� �������
#define T4    (2.35)        /// ������������ ������� ��������� ��������� 3� �������

#define LEO                 /// ��������� 6000 ������ �� ������ �����������
#define GEO                 /// ������� ����� �� ���������������


const double mu1 = 300,   mu2 = 200,   mu3 = 50;    /// ������ ������� � �������
const double u1  = 4000,  u2  = 3000,  u3  = 1750;  /// �������� ��������� ����� (� ������� ���� u3 = 115 �/� )
const double m1  = 80000, m2  = 30000, m3  = 2000;  /// ����� ����� �������
const double mp1 = 30000, mp2 = 2000,  mp3 = 5;     /// �������� ��������
      double mt1,         mt2,         mt3;         /// ����� ����
      double mf1,         mf2,         mf3;         /// ����� �������
      double thrust1,     thrust2,     thrust3;     /// ���� ��������� �������

double rho ( double h ){
/// ��������� ��������� �� ������
#define RHO0 (1.2)
#define H0 (9000.0)
#define H1 (10*H0)
    if ( h >= H1 ) return 0.;
    else return RHO0 * exp ( -h / H0 );
}

double drag ( double h, double v ){
/// ���� ������������� ������� �� ������ � ��������
#define Cx (0.4)
#define S  (2.0)
    return Cx * S * rho(h) * v*v / 2.0;
}

double mg ( double h, double m ){
/// ���� ������� �� ������ � ����� ������
    double r ;
    r = Re + h ;
    return G * Me * m / (r*r) ;
}

int i ;         /// ���������� ����� ������� � ������ �����
double  m,      /// ������� ����� ������
        ax, ay, /// ������� ��������� �� ����
        vx, vy, /// ������� �������� �� ����
        x, y,   /// ������� ���������� ������ �� ����
        v, h,   /// ��������, ������
        t,      /// ����� � ������ �����
        r, phi, /// ������ � ���� (� �������� �����������)
        mgx, mgy, /// ���� ������� �� ����
        alpha,  /// ���� �������� ������� �������� �������� ������
        xi;     /// ���� ����� �������� �������� � ������������ �� ����� ����� (� ��������)
FILE *fp , *fm ; /// ����� ��� ������ ������� ���������� � ���������

int report ( ) {
/// ������ ����������, ������������ � �������
    const char format[] =
     "%10d;%.3lf;%9.0lf;%9.0lf;% 8.2lf;%.03f;%.3lf;%.0lf;%.0lf;%.0lf\n" ;
    h = r - Re ;
    xi = acos ( (x*vx + y*vy ) / ( r * v) ) *180./M_PI ;
/// xi - ��� ���� ����� �������� �������� � ������������ �� ����� �����
/// (��������� ������������ ��������)

    fprintf ( fp, format, i, t, h, r, phi*180/M_PI, v, xi, m, x, y );
    return 0 ;
}

int message ( char *text ) {
/// ������ ���������, ������������ � ���� ���������
    fprintf ( fm , "%12.2lf:\t%s" , t , text ) ;
    return 0 ;
}

int main ( ) {
/// ����� �����:
    mt1 = K * m1;
    mt2 = K * m2;
    mt3 = K * m3;
/// ����� �������:
    mf1 = m1 - mp1 - mt1;
    mf2 = m2 - mp2 - mt2;
    mf3 = m3 - mp3 - mt3;
/// ���� ����������:
    thrust1 = u1 * mu1;
    thrust2 = u2 * mu2;
    thrust3 = u3 * mu3;

    m = m1;                 /// ��������� ����� ������
    vx = 0.;                /// X - ��� ������������
    vy = Omegae * Re ;      /// Y - ��� ����� ��������� � ����� ������
                            /// ��������� �������� �� Y �������� ��������� �����
    h = 0.;                 /// ��������� ������ = 0
    x = Re + h ;            /// ��������� ���������� - ����� �����
    y = 0.;
    t = 0.;                 /// ����� ����� � ��������

    fp = fopen ( "1.csv" , "w" ) ;
    assert ( fp ) ;
    fprintf ( fp , "i;t;h;r;phi;v;xi;m;x;y\n") ;
    fm = fopen ( "1.txt" , "w" ) ;
    assert ( fm ) ;
    int ticker = TICKER ;

    message ( "START. Stage 1 engine ON\n" ) ;

/* ---  ������������ ���� �� 1-� ������� --- */
    for ( i = 0 ; mf1 > 0; i++ ){
/// �������� �� ������ ��������� ������� 1-� ��������
        r = sqrt ( x*x + y*y ) ;
        phi = atan2 ( y, x ) ;
        h = x - Re;
        mgx = mg ( h, m ) * cos (phi) ;
        mgy = mg ( h, m ) * sin (phi) ;

/// � ��������� �� � ��������� ���� ���������, �������������
/// ������� � �������� ���� ������� �� ��� X
        ax = (thrust1 - drag ( h, v ) - mgx ) / m ;

/// � ��������� �� Y ��������� ������ �������� mg �� ��� �
        ay = -mgy / m ;

        vx += ax * DT ;     /// ��������� ��������
        vy += ay * DT ;
        v = sqrt ( vy*vy + vx*vx ) ;
        x += vx * DT ;      /// ��������� ����������
        y += vy * DT ;
        t += DT ;           /// ��������� �����
        m -= mu1 * DT ;     /// ��������� �����
        mf1 -= mu1 * DT ;
        if ( i % ticker == 0 ) report ( ) ;     /// ����� ������ � �������
    }
    message ( "Stage 1 fuel tank empty\n" ) ;
    m = m2 ;        /// ���������� ��� 1-� �������
    message ( "Stage 1 fuel tank dropped off\n") ;
    report ( ) ;
/* ---  ������������ ���� �� 1-� ������� �������� --- */

/* --- ����������� ���� ������������ ������ �������� �� ���� alpha --- */
    alpha = ALPHA;
          {
            char mes[100] ;
            sprintf ( mes , "Speed vector turn on ALHPA %.2lf deg\n" , alpha ) ;
            message ( mes ) ;
          }
    alpha *= M_PI / 180. ;
    v = sqrt ( vy*vy + vx*vx );
    alpha += atan2 ( vy , vx )  ;
    vx = v * cos ( alpha ) ;
    vy = v * sin ( alpha ) ;
    report ( ) ;

/* --- ��������� ��������� 2� ������� -- */
    message ( "Stage 2 engine ON\n" ) ;
    for ( ; mf2 > 0. ; i++ ) {   /// �� ������ ��������� �������
        double v_new ;

        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust2 / m ) * DT ;
        vx *= v_new / v ;        /// ��������� ���������� 2 �������
        vy *= v_new / v ;
        v = v_new ;

        m -= mu2 * DT ;         /// ���������� ����� ������ � �������
        mf2 -= mu2 * DT ;       /// � ���� 2� ������� �� ���� ������ ���������

        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������

        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
    }

    report ( ) ;
    message ( "Stage 2 engine fuel tank empty\n") ;
    m = m3 ;
    message ( "Stage 2 fuel tank dropped off\n") ;
    report ( ) ;

    ticker *= 100 ;      /// ����� �������� � ������� �� ���� ������� �����

#ifdef LEO
/* --- ���� �� ������ ����������� ������ -- */
    message ( "LEO flight begin\n") ;

    for ( double timer = t + 6000 ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����������
        phi = atan2 ( y, x );

        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������
        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) {
                                report ( ) ;
                                if ( t > timer ) break ;
        }
        if ( r < Re + 100000 ) message ( "Crash\n" ) ;
    }
    message ( "LEO flight finish\n") ;
    report ( ) ;
#endif

/* --- ������ ��������� ��������� 3� ������� -- */
    m = m3 ;
          {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine ON first time on %.3lf sec\n" , T3 ) ;
            message ( mes ) ;
          }
    report ( ) ;

    for ( double timer = t + T3 ;; i++ ) {
        double v_new ;

        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;        /// ��������� ���������� 3 �������
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;         /// ���������� ����� ������ � �������
        mf3 -= mu3 * DT ;       /// � ���� 3� ������� �� ���� ������ ���������

        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������

        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
        if ( t >= timer ) break ;
    }

    message ( "Stage 3 engine OFF\n") ;
    report ( ) ;


/* --- ���� � ����������� ���������� 3� ������� -- */
    message ( "Free flight begin\n") ;

    for ( ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����������
        phi = atan2 ( y, x );
        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������
        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        xi = acos ( (x*vx + y*vy ) / ( r * v) ) * 180. / M_PI ;
        if ( i % ticker == 0 ) report ( ) ;

        if ( xi >= 90. ) {
    /// ��������� ��������� ���� ����� �������� ������ ������������� ������
            message ( "Ellyspe apogee, xi = 90 deg\n" ) ;
            report ( ) ;
            break ;
        }
    }

/* --- ������ ��������� ��������� 3� ������� -- */
          {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine ON second time on %.3lf sec\n" , T4 ) ;
            message ( mes ) ;
          }
    report ( ) ;

    for ( double timer = t + T4 ; ; i++ ) {
        double v_new ;

        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;        /// ��������� ���������� 3 �������
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;         /// ���������� ����� ������ � �������
        mf3 -= mu3 * DT ;       /// � ���� 3� ������� �� ���� ������ ���������

        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������

        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
        if ( t >= timer || mf3 <= 0. ) break ;
    }

        {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine OFF. Fuel in tank %.3lf kg\n" , mf3 ) ;
            message ( mes ) ;
          }

    message ( "GOAL\n") ;
    report ( ) ;

    if ( mf3 <= 0. ) {   /// ���� ������� � 3 ���� ���������
        message ( "Stage 3 engine fuel tank empty\n") ;
        m = mp3 ;
        message ( "Stage 3 fuel tank dropped off\n") ;
        report ( ) ;
    }

    ticker *= 2 ;

#ifdef GEO
/* --- ���� � ����������� ���������� 3� ������� �� ��������������� -- */
    message ( "Geostat flight begin\n") ;

    for ( double timer = t + Day * 1.05 ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// ��������� ������ � �������� �����������
        phi = atan2 ( y, x );

        mgx = mg ( h, m ) * cos (phi); /// �������� ���� ������� �� ���
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// �������� ��������� �� ���
        ay = -mgy / m; /// ������������� ������� ��� �� ���������
        vx += ax * DT ;  /// �������� ���������, ��������� � �������
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) {
                                report ( ) ;
                                if ( t > timer )  break ;
                                }
    }

    message ( "MISSION COMPLETED\n"); /// ��������� ���������
    report ( ) ;

#endif

 fclose ( fm ) ;
 fclose ( fp ) ;
 return 0 ;
}

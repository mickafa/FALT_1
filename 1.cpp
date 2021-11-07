#include <bits/stdc++.h>
using namespace std ;
#define Day (86164.09054)   /// длительность суток
#define Me (5.9726e24)      /// масса Земли
#define Re (6378100.0)      /// экваториальный радиус
#define G (6.6743015e-11)   /// гравитационная постоянная
#define DT (0.001)          /// шаг по времени
#define TICKER (10000)    /// начальная периодичность репорта параметров в таблицу (в тиках)

#define Omegae (2*M_PI/Day) /// угловая скорость Земли
#define Ro (42165188.26463) /// радиус геостационарной орбиты
#define Vo (Omegae*Ro)      /// орбитальная скорость на геостационарной

#define K (0.05)            /// удельная масса баков

#define ALPHA (59.0)        /// угол поворота вектора скорости внешней силой
#define T3    (35.335)      /// длительность первого включения двигателя 3й ступени
#define T4    (2.35)        /// длительность второго включения двигателя 3й ступени

#define LEO                 /// пролететь 6000 секунд по низкой околоземной
#define GEO                 /// сделать виток по геостационарной


const double mu1 = 300,   mu2 = 200,   mu3 = 50;    /// расход топлива в секунду
const double u1  = 4000,  u2  = 3000,  u3  = 1750;  /// скорость истечения газов (в задании было u3 = 115 м/с )
const double m1  = 80000, m2  = 30000, m3  = 2000;  /// общая масса ступени
const double mp1 = 30000, mp2 = 2000,  mp3 = 5;     /// полезная нагрузка
      double mt1,         mt2,         mt3;         /// масса бака
      double mf1,         mf2,         mf3;         /// масса топлива
      double thrust1,     thrust2,     thrust3;     /// тяга двигателя ступени

double rho ( double h ){
/// плотность атмосферы от высоты
#define RHO0 (1.2)
#define H0 (9000.0)
#define H1 (10*H0)
    if ( h >= H1 ) return 0.;
    else return RHO0 * exp ( -h / H0 );
}

double drag ( double h, double v ){
/// сила сопротивления воздуха от высоты и скорости
#define Cx (0.4)
#define S  (2.0)
    return Cx * S * rho(h) * v*v / 2.0;
}

double mg ( double h, double m ){
/// сила тяжести от высоты и массы ракеты
    double r ;
    r = Re + h ;
    return G * Me * m / (r*r) ;
}

int i ;         /// количество тиков времени с начала полёта
double  m,      /// текущая масса ракеты
        ax, ay, /// текущие ускорения по осям
        vx, vy, /// текущие скорости по осям
        x, y,   /// текущие координаты ракеты по осям
        v, h,   /// скорость, высота
        t,      /// время с начала полёта
        r, phi, /// радиус и угол (в полярных координатах)
        mgx, mgy, /// сила тяжести по осям
        alpha,  /// угол поворота вектора скорости внешними силами
        xi;     /// угол между вектором скорости и направлением на центр Земли (в градусах)
FILE *fp , *fm ; /// файлы для записи таблицы параметров и сообщений

int report ( ) {
/// строка параметров, записываемая в таблицу
    const char format[] =
     "%10d;%.3lf;%9.0lf;%9.0lf;% 8.2lf;%.03f;%.3lf;%.0lf;%.0lf;%.0lf\n" ;
    h = r - Re ;
    xi = acos ( (x*vx + y*vy ) / ( r * v) ) *180./M_PI ;
/// xi - это угол между вектором скорости и направлением на центр Земли
/// (скалярное произведение векторов)

    fprintf ( fp, format, i, t, h, r, phi*180/M_PI, v, xi, m, x, y );
    return 0 ;
}

int message ( char *text ) {
/// строка сообщения, записываемая в файл сообщений
    fprintf ( fm , "%12.2lf:\t%s" , t , text ) ;
    return 0 ;
}

int main ( ) {
/// массы баков:
    mt1 = K * m1;
    mt2 = K * m2;
    mt3 = K * m3;
/// массы топлива:
    mf1 = m1 - mp1 - mt1;
    mf2 = m2 - mp2 - mt2;
    mf3 = m3 - mp3 - mt3;
/// тяга двигателей:
    thrust1 = u1 * mu1;
    thrust2 = u2 * mu2;
    thrust3 = u3 * mu3;

    m = m1;                 /// стартовая масса ракеты
    vx = 0.;                /// X - ось вертикальная
    vy = Omegae * Re ;      /// Y - ось вдоль горизонта в точке старта
                            /// начальная скорость по Y создаётся вращением Земли
    h = 0.;                 /// начальная высота = 0
    x = Re + h ;            /// начальные координаты - центр Земли
    y = 0.;
    t = 0.;                 /// время полёта в секундах

    fp = fopen ( "1.csv" , "w" ) ;
    assert ( fp ) ;
    fprintf ( fp , "i;t;h;r;phi;v;xi;m;x;y\n") ;
    fm = fopen ( "1.txt" , "w" ) ;
    assert ( fm ) ;
    int ticker = TICKER ;

    message ( "START. Stage 1 engine ON\n" ) ;

/* ---  вертикальный взлёт на 1-й ступени --- */
    for ( i = 0 ; mf1 > 0; i++ ){
/// работаем до полной выработки топлива 1-й ступенью
        r = sqrt ( x*x + y*y ) ;
        phi = atan2 ( y, x ) ;
        h = x - Re;
        mgx = mg ( h, m ) * cos (phi) ;
        mgy = mg ( h, m ) * sin (phi) ;

/// в ускорении по Х учитываем тягу двигателя, сопротивление
/// воздуха и проекцию силы тяжести на ось X
        ax = (thrust1 - drag ( h, v ) - mgx ) / m ;

/// в ускорении по Y учитываем только проекцию mg на ось Х
        ay = -mgy / m ;

        vx += ax * DT ;     /// обновляем скорости
        vy += ay * DT ;
        v = sqrt ( vy*vy + vx*vx ) ;
        x += vx * DT ;      /// обновляем координаты
        y += vy * DT ;
        t += DT ;           /// обновляем время
        m -= mu1 * DT ;     /// обновляем массы
        mf1 -= mu1 * DT ;
        if ( i % ticker == 0 ) report ( ) ;     /// пишем репорт в таблицу
    }
    message ( "Stage 1 fuel tank empty\n" ) ;
    m = m2 ;        /// сбрасываем бак 1-й ступени
    message ( "Stage 1 fuel tank dropped off\n") ;
    report ( ) ;
/* ---  вертикальный взлёт на 1-й ступени завершён --- */

/* --- посторонняя сила поворачивает вектор скорости на угол alpha --- */
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

/* --- включение двигателя 2й ступени -- */
    message ( "Stage 2 engine ON\n" ) ;
    for ( ; mf2 > 0. ; i++ ) {   /// до полной выработки топлива
        double v_new ;

        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных коорд
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust2 / m ) * DT ;
        vx *= v_new / v ;        /// ускорение двигателем 2 ступени
        vy *= v_new / v ;
        v = v_new ;

        m -= mu2 * DT ;         /// уменьшение массы ракеты и топлива
        mf2 -= mu2 * DT ;       /// в баке 2й ступени за счёт работы двигателя

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем

        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
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

    ticker *= 100 ;      /// чтобы репортов в таблице не было слишком много

#ifdef LEO
/* --- полёт по низкой околоземной орбите -- */
    message ( "LEO flight begin\n") ;

    for ( double timer = t + 6000 ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных координатах
        phi = atan2 ( y, x );

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем
        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
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

/* --- первое включение двигателя 3й ступени -- */
    m = m3 ;
          {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine ON first time on %.3lf sec\n" , T3 ) ;
            message ( mes ) ;
          }
    report ( ) ;

    for ( double timer = t + T3 ;; i++ ) {
        double v_new ;

        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных коорд
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;        /// ускорение двигателем 3 ступени
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;         /// уменьшение массы ракеты и топлива
        mf3 -= mu3 * DT ;       /// в баке 3й ступени за счёт работы двигателя

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем

        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
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


/* --- полёт с выключенным двигателем 3й ступени -- */
    message ( "Free flight begin\n") ;

    for ( ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных координатах
        phi = atan2 ( y, x );
        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем
        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        xi = acos ( (x*vx + y*vy ) / ( r * v) ) * 180. / M_PI ;
        if ( i % ticker == 0 ) report ( ) ;

        if ( xi >= 90. ) {
    /// прерываем свободный полёт когда достигли апогея эллиптической орбиты
            message ( "Ellyspe apogee, xi = 90 deg\n" ) ;
            report ( ) ;
            break ;
        }
    }

/* --- второе включение двигателя 3й ступени -- */
          {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine ON second time on %.3lf sec\n" , T4 ) ;
            message ( mes ) ;
          }
    report ( ) ;

    for ( double timer = t + T4 ; ; i++ ) {
        double v_new ;

        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных коорд
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;        /// ускорение двигателем 3 ступени
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;         /// уменьшение массы ракеты и топлива
        mf3 -= mu3 * DT ;       /// в баке 3й ступени за счёт работы двигателя

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем

        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
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

    if ( mf3 <= 0. ) {   /// если топливо в 3 баке кончилось
        message ( "Stage 3 engine fuel tank empty\n") ;
        m = mp3 ;
        message ( "Stage 3 fuel tank dropped off\n") ;
        report ( ) ;
    }

    ticker *= 2 ;

#ifdef GEO
/* --- полёт с выключенным двигателем 3й ступени по геостационарной -- */
    message ( "Geostat flight begin\n") ;

    for ( double timer = t + Day * 1.05 ; ; i++ ) {
        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных координатах
        phi = atan2 ( y, x );

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m; /// проекции ускорения на оси
        ay = -mgy / m; /// сопротивление воздуха уже не учитываем
        vx += ax * DT ;  /// пересчёт скоростей, координат и времени
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

    message ( "MISSION COMPLETED\n"); /// завершаем программу
    report ( ) ;

#endif

 fclose ( fm ) ;
 fclose ( fp ) ;
 return 0 ;
}

#include <bits/stdc++.h>
#define Day (86164.09054)   /// длительность суток (собственное вращение Земли)
#define Me (5.9726e24)      /// масса Земли
#define Re (6378100.0)      /// экваториальный радиус
#define G (6.6743015e-11)   /// гравитационная постоянная
#define DT (0.001)          /// шаг расчета по времени
#define TICKER (10000)      /// периодичность репорта параметров в таблицу (в тиках)

#define Omegae (2*M_PI/Day) /// угловая скорость Земли
#define Ro (42165188.26463) /// радиус геостационарной орбиты
#define Vo (Omegae*Ro)      /// орбитальная скорость на геостационарной

#define K (0.05)            /// удельная масса баков

#define ALPHA (60.136)      /// угол поворота вектора скорости внешней силой
#define T3    (35.312)      /// длительность первого включения двигателя 3й ступени

/// #define LEO              /// сделать виток по низкой околоземной - НЕ ДЕЛАЕМ
#define GEO                  /// сделать виток по геостационарной

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
        xi;     /// угол между вектором скорости и отвесом от центра Земли (в градусах)
FILE *fp , *fm ; /// файлы для записи таблицы параметров и сообщений

int report ( ) {
/// репорт - строка параметров, записываемая в таблицу
    static char format[] =
     "%10d;%.3lf;%9.0lf;%9.0lf;% 8.2lf;%.03f;%.3lf;%.0lf;%.0lf;%.0lf\n" ;
    h = r - Re ;
    xi = acos ( (x*vx + y*vy ) / ( r * v) ) * 180./M_PI ;
/// xi - это угол между вектором скорости и отвесом от центра Земли
/// (через скалярное произведение векторов)

    fprintf ( fp, format, i, t, h, r, phi*180/M_PI, v, xi, m, x, y );
     printf (     format, i, t, h, r, phi*180/M_PI, v, xi, m, x, y );
    return 0 ;
}

int message ( const char *text ) {
/// строка сообщения, записываемая в файл сообщений
    fprintf ( fm , "%12.2lf:\t%s" , t , text ) ;
     printf (      "%12.2lf:\t%s" , t , text ) ;
    return 0 ;
}

int main ( ) {
/// устанавливаем начальные условия:
/// массы баков
    mt1 = K * m1;
    mt2 = K * m2;
    mt3 = K * m3;
/// массы топлива
    mf1 = m1 - mp1 - mt1;
    mf2 = m2 - mp2 - mt2;
    mf3 = m3 - mp3 - mt3;
/// тяга двигателей
    thrust1 = u1 * mu1;
    thrust2 = u2 * mu2;
    thrust3 = u3 * mu3;

    m = m1;                 /// стартовая масса ракеты
    vx = 0.;                /// X - ось вертикальная
    vy = Omegae * Re ;      /// Y - ось вдоль горизонта в точке старта
                            /// начальная скорость по Y создаётся вращением Земли
    h = 0.;                 /// начальная высота = 0
    x = Re + h ;            /// начальные координаты в системе центра Земли
    y = 0.;
    t = 0.;                 /// время полёта в секундах


    fp = fopen ( "1.csv" , "w" ) ; /// открываем файлы для записи репортов,
    assert ( fp ) ;
    fprintf ( fp , "i;t;h;r;phi;v;xi;m;x;y\n") ;
    fm = fopen ( "1.txt" , "w" ) ; /// и сообщений
    assert ( fm ) ;
    fprintf ( fm , "Flight time:\tMessage:\n") ;
    int ticker = TICKER ;           /// периодичность записи репортов: 10 секунд

    message ( "START. Stage 1 engine ON\n" ) ;   /// Поехали

/* ---  вертикальный взлёт на 1-й ступени --- */
    for ( i = 0 ; mf1 > 0; i++ ){
/// работаем до полной выработки топлива mf1 из бака 1-й ступени
/// пересчитываем текущие координаты в полярные
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
        m -= mu1 * DT ;     /// уменьшаем массу ракеты
        mf1 -= mu1 * DT ;   /// и топлива в баке 1-й ступени
        if ( i % ticker == 0 ) report ( ) ;     /// пишем репорт в таблицу
    }
    message ( "Stage 1 fuel tank empty\n" ) ;
    m = m2 ;                /// сбрасываем бак 1-й ступени
    message ( "Stage 1 fuel tank dropped off\n") ;
    report ( ) ;
/* ---  вертикальный взлёт на 1-й ступени завершён --- */

/* --- высшая сила поворачивает вектор скорости на угол alpha --- */
    alpha = ALPHA;
          {
            char mes[100] ;
            sprintf ( mes , "Speed vector turn on ALHPA %.3lf deg\n" , alpha ) ;
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
    for ( ; mf2 > 0. ; i++ ) {   /// до полной выработки топлива mf2
        double v_new ;

        r = sqrt ( x*x + y*y );  /// положение ракеты в полярных коорд
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust2 / m ) * DT ;
        vx *= v_new / v ;        /// удлинняем вектор скорости в соответствии с тягой двигателя
        vy *= v_new / v ;        /// и пропорционально распределяем новую скорость по осям
        v = v_new ;              /// скорость, пересчитанная после добавления импульса тяги

        m -= mu2 * DT ;         /// уменьшение массы ракеты
        mf2 -= mu2 * DT ;       /// и топлива в баке 2й ступени

        mgx = mg ( h, m ) * cos (phi); /// проекции силы тяжести на оси
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m;          /// проекции ускорения от силы тяжести на оси
        ay = -mgy / m;          /// сопротивление воздуха уже не учитываем

        vx += ax * DT ;         /// пересчёт скоростей за счёт действия силы тяжести
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;          /// обновление координат и времени
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
    }
    report ( ) ;
    message ( "Stage 2 engine fuel tank empty\n") ;
    m = m3 ;
    message ( "Stage 2 fuel tank dropped off\n") ;
    report ( ) ;
/* ---  вторая ступень отработала --- */


    ticker *= 100 ;     /// периодичность увеличиваем в 100 раз до 1000 с,
                        /// чтобы не забивать таблицу почти одинаковыми репортами

#ifdef LEO
/* --- полёт по низкой околоземной орбите, ОТМЕНЁН -- */
    message ( "LEO flight begin\n") ;
cout <<  (2.*M_PI*r/v) * 1. << '*' << t << endl ;
    for ( double timer = t + (2.*M_PI*r/v) * 1. ; t < timer ; i++ ) {

        r = sqrt ( x*x + y*y );  /// летим по таймеру 1 виток
        phi = atan2 ( y, x );

        mgx = mg ( h, m ) * cos (phi);
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m;
        ay = -mgy / m;
        vx += ax * DT ;
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
        if ( r < Re + 100000 ) message ( "Crash\n" ) ;
    }
    message ( "LEO flight finish\n") ;
    report ( ) ;
/* --- конец полёта по низкой околоземной орбите -- */
#endif


/* --- первое включение двигателя 3й ступени, уход с низкой околоземной -- */
          {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine ON first time on %.3lf sec\n" , T3 ) ;
            message ( mes ) ;
          }
    report ( ) ;
    for ( double timer = t + T3 ; t < timer ; i++ ) {  /// по таймеру на время T3
        double v_new ;

        r = sqrt ( x*x + y*y );  /// алгоритм аналогичен как был для 2й ступени
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;        /// ускорение двигателем 3 ступени
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;
        mf3 -= mu3 * DT ;

        mgx = mg ( h, m ) * cos (phi);
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m;
        ay = -mgy / m;

        vx += ax * DT ;
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
    }

        {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine OFF. Fuel in tank %.3lf kg\n" , mf3 ) ;
            message ( mes ) ;
          }    report ( ) ;
/* --- конец первого включения двигателя 3й ступени -- */


/* --- полёт с выключенным двигателем по переходной эллиптической орбите -- */
    message ( "Free flight begin\n") ;
    for ( ;  xi < 90.0 ; i++ ) {  /// до прибытия в апогей: пока угол xi не станет 90
        r = sqrt ( x*x + y*y );
        phi = atan2 ( y, x );
        mgx = mg ( h, m ) * cos (phi);
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m;
        ay = -mgy / m;
        vx += ax * DT ;
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        xi = acos ( (x*vx + y*vy ) / ( r * v) ) * 180. / M_PI ;
        if ( i % ticker == 0 ) report ( ) ;
    }
    message ( "Ellyspe apogee, xi = 90 deg\n" ) ;
    report ( ) ;
/* --- мы в апогее переходной орбиты -- */


/* --- второе включение двигателя 3й ступени для стабилизации на геостационарной -- */
    message ( "Stage 3 engine ON second time\n" ) ;
    report ( ) ;
    for (  ; v < Vo ; i++ ) {  /// работаем пока не достигнута орбитальная скорость
        double v_new ;

        r = sqrt ( x*x + y*y );
        phi = atan2 ( y, x );
        v_new = sqrt ( vy*vy + vx*vx ) + ( thrust3 / m ) * DT ;
        vx *= v_new / v ;
        vy *= v_new / v ;
        v = v_new ;

        m -= mu3 * DT ;
        mf3 -= mu3 * DT ;

        mgx = mg ( h, m ) * cos (phi);
        mgy = mg ( h, m ) * sin (phi);

        ax = -mgx / m;
        ay = -mgy / m;

        vx += ax * DT ;
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
    }
        {
            char mes[100] ;
            sprintf ( mes , "Stage 3 engine OFF. Fuel in tank %.3lf kg\n" , mf3 ) ;
            message ( mes ) ;
          }
    report ( ) ;
/* --- собственно, мы на месте -- */

#ifdef GEO
/* --- контрольный виток по геостационарной -- */
    ticker *= 2 ;       /// увеличиваем периодичность репортов в 2 раза до 2000 секунд
    message ( "GEO flight begin\n") ;

    int flag = 0 ;
    for( double timer = t + Day ; t < timer ; i++ ) {///по таймеру на время одного витка
        r = sqrt ( x*x + y*y );
        phi = atan2 ( y, x );
        mgx = mg ( h, m ) * cos (phi);
        mgy = mg ( h, m ) * sin (phi);
        ax = -mgx / m;
        ay = -mgy / m;
        vx += ax * DT ;
        vy += ay * DT ;
        v = sqrt ( vx*vx + vy*vy ) ;
        x += vx * DT ;
        y += vy * DT ;
          
        if ( y >= 0. && flag == 0 ) {
            char mes[100] ;
            sprintf ( mes , "Zero point: x=%.0lf; v=%.0lf\n" , x, v ) ;
            message ( mes ) ;
            report ( ) ;
            flag = 1 ;
        }

        t += DT ;
        if ( i % ticker == 0 ) report ( ) ;
    }

    message ( "GEO flight completed\n"); /// завершаем программу
    report ( ) ;
/* --- контрольный виток по геостационарной завершён -- */
#endif

        {
            char mes[100] ;
            sprintf ( mes , "GOAL orbit:\tRo= %.0lf m\tVo=%.2lf\n" , Ro, Vo ) ;
            message ( mes ) ;
            sprintf ( mes , "Reached orbit:\tr = %.0lf m\t v=%.2lf\n" , r, v ) ;
            message ( mes ) ;

          }

 fclose ( fm ) ;
 fclose ( fp ) ;
 return 0 ;
}

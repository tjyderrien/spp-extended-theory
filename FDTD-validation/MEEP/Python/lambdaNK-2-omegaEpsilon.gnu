#!gnuplot
c=3E8
pi=acos(-1)

epsilon(n,k)=(n*{1,0}+k*{0,1})**2
omega(lambda)=2*pi*c/lambda
xscale=1E-6

set xlabel 'Frequency (s^{-1})'
set ylabel 'Epsilon'

set xrange [1E9:1E17]

set log xy

plot "AgRakic.dat" u (omega($1*xscale)):(real(epsilon($2,$3))) w l t 'Re(+epsilon)', \
"AgRakic.dat" u (omega($1*xscale)):(-real(epsilon($2,$3))) w l t 'Re(-epsilon)', \
"AgRakic.dat" u (omega($1*xscale)):(imag(epsilon($2,$3))) w l t 'Im(epsilon)'

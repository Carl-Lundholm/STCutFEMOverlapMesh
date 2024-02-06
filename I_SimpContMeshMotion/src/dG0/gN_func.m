function g = gN_func(x,t)

global x0_init x0_fin

g_x0init = 0;
g_x0fin = 0;

g_x0init = g_x0init + (t - t);
g_x0fin = g_x0fin + (t - t);

g = 1/(x0_fin - x0_init)*((x0_fin - x)*g_x0init + (x - x0_init)*g_x0fin);


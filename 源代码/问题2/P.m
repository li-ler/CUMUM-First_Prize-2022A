function ans = P(c, n, delta_t)
    w = 2.2143;
    m_a1 = 1165.992;
    c1 = 167.8395;
    f = 4890;

    m_f = 4866;
    r = 1;
    m_o = 2433;
    rho = 1025;
    k = 80000;
    g = 9.8;
    l0 = 0.5;

    xf0 = -2;
    xo0 = -1.8;

    %c = 10000;

    T = 2 * pi / w;

    dp_c=@(t,p,c)[p(2); 
              (-(c1 + c) * p(2) - k * p(1) + c * p(4) + k * p(3) - k * l0 + ...
              rho * g * V(p(1)) - m_f * g +...
              f * cos(w * t)) / (m_a1 + m_f);
              p(4); 
              (c * p(4) + k * p(3) - c * p(2) - k * p(1) + m_o * g - k * l0) / (-m_o);];
    dp1 = @(t,p)dp_c(t, p, c);
    dp2 = @(t,p)dp_c(t, p, c*((abs(p(2) - p(4)))^n));
    sol1=ode45(dp2,[0,40*T],[-2 0 -1.8 0]);
    t = 0:delta_t:40*T;
    p1 = deval(sol1,t);
    xf1 = p1(1,:);
    vf1 = p1(2,:);
    xo1 = p1(3,:);
    vo1 = p1(4,:);
    p =  c * (vo1 - vf1) .* (vo1 - vf1);
    ans = sum(p) * delta_t / (40*T);
end
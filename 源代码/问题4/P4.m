function ans = P(c, ct, delta_t)
    w = 1.9806;
    m_a1 = 1091.099;
    Ia = 7142.493;
    c1 = 528.5018;
    cp = 1655.909;
    f = 1760;
    L = 2140;

    m_f = 4866;
    Rf = 1;
    Hf1 = 3;
    Hf2 = 0.8;
    Sf1 = pi*Rf^2;
    Sf2 = pi*Rf^2;
    Sf3 = 2*pi*Rf*Hf1;
    l = sqrt(Hf2^2 + Rf^2);
    Sf4 = pi*Rf*l;
    mf1 = Sf1 / (Sf1 + Sf2 + Sf3 + Sf4) * m_f;
    mf2 = Sf2 / (Sf1 + Sf2 + Sf3 + Sf4) * m_f;
    mf3 = Sf3 / (Sf1 + Sf2 + Sf3 + Sf4) * m_f;
    mf4 = Sf4 / (Sf1 + Sf2 + Sf3 + Sf4) * m_f;
    If1 = 1 / 4 * mf1 * Rf^2 + mf1 * Hf1^2;
    If2 = 1 / 4 * mf2 * Rf^2;
    If3 = 1 / 2 * mf3 * Rf^2 + 1 / 12 * mf3 * Hf1^2 + mf3 * (Hf1 / 2)^2;
    dh = [0.01:0.01:Hf2];
    dr = dh * Rf / Hf2;
    dl = 2 * pi * dr;
    dm = mf4 / (sum(dl)) * dl;
    If4 = 1 / 2 * dm.*(dr.^2) + dm.*((Hf2 - dh).^2);
    If4 = sum(If4);

    If = If1 + If2 + If3 + If4;

    m_o = 2433;
    Ro = 0.5;
    Ho = 0.5;
    rho = 1025;
    ks = 80000;
    g = 9.8;
    l0 = 0.5;
    kt = 250000;
    Khs = 8890.7;

    xf0 = -2;
    xo0 = -1.8;


    T = 2 * pi / w;

    dp_c=@(t,p,c)[p(2); 
              (f*cos(w*t) - c1*p(2) + (rho*g*V(p(1)) - m_f*g) * cos(p(5)) +...
              (-ks*(l0 - (p(3) - p(1))) + c * (p(4) - p(2))) * cos(p(7)-p(5))) / (m_f + m_a1);
              p(4); 
              (ks*(l0 - (p(3)-p(1))) - c * (p(4)-p(2)) - m_o*g*cos(p(7))) / m_o;
              p(6);
              (L * cos(w*t) - cp*p(6) - Khs*p(5) + kt*(p(7)-p(5)) + ct*(p(8)-p(6)) -...
              m_o*g*sin(p(7))*(Ho / 2 + p(3)-p(1))) / (If + Ia);
              p(8);
              (-kt*(p(7)-p(5)) - ct*(p(8)-p(6)) + m_o*g*sin(p(7))*(Ho / 2 + p(3)-p(1))) /...
              (1 / 12 * m_o * (3*Ro^2 + Ho^2) + m_o*(Ho/2 + p(3) - p(1))^2)];
    dp1 = @(t,p)dp_c(t, p, c);

    sol1=ode45(dp1,[0,40*T],[-2 0 -1.8 0 0 0 0 0]);
    t = 0:delta_t:40*T;
    p1 = deval(sol1,t);
    xf1 = p1(1,:);
    vf1 = p1(2,:);
    xo1 = p1(3,:);
    vo1 = p1(4,:);
    theta_f1 = p1(5,:);
    wf1 = p1(6,:);
    theta_o1 = p1(7,:);
    wo1 = p1(8,:);
    p =  c * (vo1 - vf1) .* (vo1 - vf1) + ct * (wo1 - wf1) .* (wo1 - wf1);
    ans = sum(p) * delta_t / (40*T);
end
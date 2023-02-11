function Im = generate_regular_polyhedron(s, c, r, N)

    % s = [512,512]; %size
    % c = s/2; %center
    % r = 200; %radius
    % %approximate circle with a regular polygon
    % N = 100;

    xs = zeros(1,N);
    ys = zeros(1,N);
    xs(1) = c(2);
    ys(1) = c(1)-r;
    for i = 2:N
        xs(i) = xs(i-1)+2*r*sin(pi/N)*cos(pi/N+(i-2)*2*pi/N);
        ys(i) = ys(i-1)+2*r*sin(pi/N)*sin(pi/N+(i-2)*2*pi/N);
    end

    Im = poly2mask(xs,ys,s(1),s(2));

end
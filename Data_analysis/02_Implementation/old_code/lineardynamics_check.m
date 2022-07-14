start_pos = [-0.1, 0.4]; 
% start_pos = [-0.045, 0.14]; 
target = [0.4, 0.1];
current = start_pos; 
dt = 0.0742;
t = 0 ;
trace = [current];
while t < 10
    v = -5 * (current- target);  % x' = -10(x-0.4)
    current = current + dt * v;
    trace = [trace; current];
    t = t+ dt;
end

scatter(trace(:,1), trace(:,2),8, 'g');
hold on
circle(0.4,0.1,0.02,2,'c');
hold on
plot(trace(:,1), trace(:,2));
hold on
plot(0.4,0.1, 'x', 'color','r', 'LineWidth',10)
grid on
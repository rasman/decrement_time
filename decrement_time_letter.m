warning('off','MATLAB:print:FigureTooLargeForPage')

mPatient = patBuilder(50,70, 170, 1);
dur = 4*60;
plan = [4, 0];
time = 0:1/60:dur;
t1 = 0:1/60:600;
u1 = zeros(size(t1));
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
set(h,'Position',[100 100 1250 500])

filename = 'decrement_time_pk_var_print.tiff';

input =  build_TCI_eleveld (mPatient, plan, dur, 0);
val = dur*60;

N = 2000;
location = nan(N,5);
d_time = nan(N,5);
start_con = nan(N,1);

varysize = 1.5;

for iter = 1:N
    sys_ele = eleveld18_vary(mPatient,varysize);
    sys_ele = sys_ele(2);
    [y1,tOut,x] = lsim(sys_ele,input,time);
    y = lsim(sys_ele,u1,t1, x(val,:));

    start_con(iter) = y1(end);
    
    d_time(iter, 1) = find(y < y1(end)*.75,1);
    d_time(iter, 2) = find(y < y1(end)*.625,1);
    d_time(iter, 3) = find(y < y1(end)*.5,1);
    d_time(iter, 4) = find(y < y1(end)*.375,1);
    d_time(iter, 5) = find(y < y1(end)*.25,1);

    % if y1(end) <= 3e-3
    try
        location(iter, 5) = find(y < 1e-3,1);
        location(iter, 4) = find(y < 1.5e-3,1);
        location(iter, 3) = find(y < 2e-3,1);
        location(iter, 2) = find(y < 2.5e-3,1);
        location(iter, 1) = find(y < 3e-3,1);      
    catch
        disp('fail');
    end
    if mod(iter,100) ==0
        disp(iter)
    end
end



subplot(1,2,1)

mega = t1(location);
boxplot(mega,{'25%','37.5%','50%','62.5%','75%'},'Colors','rkbmc','Symbol','ko')


ylim([0 50])
ylim([0 50])
xlabel('Decrement percentage')
ylabel('Time to reach decrement value (min)')
title ('Decrement time relative to set concentration')



subplot(1,2,2)

mega2 = t1(d_time);
boxplot(mega2,{'25%','37.5%','50%','62.5%','75%'},'Colors','rkbmc','Symbol','ko')

ylim([0 50])
xlabel('Decrement percentage')
ylabel('Time to reach decrement value (min)')
title ('Decrement time relative to final concentration')

frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'tiff','Resolution',300);
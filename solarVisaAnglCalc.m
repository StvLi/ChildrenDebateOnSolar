%   关于太阳视角讨论的仿真程序
%   2023-07-01 StvLi 没事闲的
clear;
clc;
%   问题描述：
%       与章QF、赵HN同学讨论以下问题：
%       太阳在正午视角大？还是在日出时视角大？
%
%   仿真模型：
%       几何数据：地球半径 6.4e3 太阳半径7.0e5 日地距离1.5e8
%       折射模型：单层模型 1.0e2 大气 1.0003 真空 1.0000 
%       模型描述：
%           在日心、地心、观测者三点平面内建立几何光学模型
%           地心位于(0,0)处 日心位于(1.5e8,0)处
%           地面为以(0,0)为圆心 6400为半径的圆
%           单层大气折射界面为以(0,0)为圆心 6500为半径的圆
%           观测者位于地面上 大气中的“视线”过观测者
%           观测者位置通过日心-地心连线与观测者-地心连线夹角描述
%           大气中的“视线”在折射界面处发生斯涅耳折射
%           真空中的“视线”距离(1.5e8,0)小于7.0e5即认为可见太阳
%           对以上 观测者位置 和 视角 依次遍历 (x轴正向)

%   数据初始化
%       模型数据
x_earthCenter = 0;
y_earthCenter = 0;
r_earth = 6.4e+3;
r_atmos = 6.5e+3;
n_atmos = 1.0003;
n_vacuo = 1.0000;
x_solarCenter = 1.5e8;
y_solarCenter = 0;
r_solar = 7.0e5;
%       求解设定
obseLoca_lower = -90 * (pi/180);
obseLoca_upper = +90 * (pi/180);
obseLoca_step  = 0.5 * (pi/180);
optiAngl_lower = -1e-2;
optiAngl_upper = +1e-2;
optiAngl_step  = 5e-7;
%       记录序号
count = 1;

%   遍历 观测者位置：0°(正午)~+90°(日出？其实不完全是)
for obseLoca = obseLoca_lower : obseLoca_step : obseLoca_upper
    disp("routine")
    disp(count)
    flag_pre = 0;
    flag_now = 0;
    visaAngl_min(count) = optiAngl_lower;
    visaAngl_max(count) = optiAngl_upper;
%       遍历 大气内光线角度：-1e-2rad ~+1e-2rad
    for optiAngl = optiAngl_lower : optiAngl_step : optiAngl_upper
%           计算入射点
        a = tan(optiAngl);
        b = - a * cos(obseLoca) * r_earth + sin(obseLoca) * r_earth;
        delta = 4 * a*a * b*b + 4 * (a*a + 1)*(r_atmos*r_atmos - b*b);
        x_point = (-2*a*b + sqrt(delta))/(2*a*a + 2);
        y_point = a * x_point + b;
%           计算折射入射角
        theta_in = - atan(y_point/x_point) + atan((y_point - ...
            sin(obseLoca) * r_earth)/(x_point - cos(obseLoca) * r_earth));
%           计算折射出射角
        theta_out = asin((n_atmos / n_vacuo) * sin(theta_in));
%           判断是否可见太阳
        A = tan (theta_out + atan(y_point/x_point));
        B = -1;
        C = -A*cos(obseLoca)*r_earth - B*sin(obseLoca)*r_earth;
        distance = abs(A*x_solarCenter + B*y_solarCenter + C) / ...
            sqrt(A^2 + B^2);
        flag_pre = flag_now;
        if distance < r_solar
            flag_now = 1;
        else
            flag_now = 0;
        end
%       记录 观测者位置-太阳视角大小
        if flag_pre == 0 && flag_now == 1
            visaAngl_min(count) = optiAngl;
        end
        if flag_pre == 1 && flag_now == 0
            visaAngl_max(count) = optiAngl;
        end
    end
    if (visaAngl_max(count) - visaAngl_min(count)) > 0
        visaAnglDegr(count) = visaAngl_max(count) - visaAngl_min(count);
    end
    count = count+1;
end
%   绘制 观测者位置-太阳视角大小 图
plot(visaAnglDegr)


%% loop
step_t=3000;
N=1001;
R_num=1;
T_ev_cleo=zeros(N,step_t+2,3);
I_ev_cleo=zeros(step_t+2,1,3);
Vc_ev_cleo=zeros(step_t+2,1,3);
R_ev_cleo=zeros(step_t+2,1,3);
for Rs=[0,150,350]%serial resistance

    T_ev=zeros(N,step_t+2);
    I_ev=zeros(step_t+2,1);
    R_ev=zeros(step_t+2,1);
    Vc_ev=zeros(step_t+2,1);
    
    % Electrical parameters
    Lk=750*10^(-9);%807.7;%nH
    Z0=50;%ohm
    Cbt=100*10^(-9);%nF
    
    %SNSPD parameters
    Tc=10;%K
    Tsub=2.5;%K
    Delta_t=1*10^(-12);%ps
    Delta_x=3.75*10^(-9);%nm
    Ic0=13.13*10^(-6);%uA
    Ibias=0.9*Ic(Tsub,Ic0,Tc);
    
    %initial state
    T_ev(:,1)=ones(N,1)*Tsub;
    T_0=ones(N,1)*Tsub;
    for j=499:503%initial hot spot size
        T_0(j)=sqrt(1-sqrt(Ibias/Ic0))*Tc;
    end
    T_ev(:,2)=T_0;
    I_ev(1:2)=Ibias;
    Vc_ev(1:2)=Ibias*(R_ev(1:2)+Rs);
    
    %geometry of the SNSPD   
    w=70*10^(-9);%width
    d=10*10^(-9);%thickness
    Rd=300;%ohm%sheet resistance
    R_ev(2)=Rd/w*15*10^(-9);%15nm is the hotspot size, corresponding to the
    % for loop T_0(j) for the initial state
    
    for i=3:step_t+2
        
        T_0=T_ev(:,i-1);
        I_0=I_ev(i-1);
        T_ini=T_ev(:,i-2);
        I_00=I_ev(i-2);
    
        A_c0=linspace(1,N,N);
        A_r0=linspace(1,N,N);
        A_0=ones(N,1);
        A_0(2:N-1)=1+Delta_t/2/d*alpha(T_0(2:N-1))./ct(T_0(2:N-1),Tc,I_0,Ic0)+...
            2*Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        A_c_left=linspace(2,N-1,N-2);
        A_r_left=linspace(2,N-1,N-2)-1;
        A_0_left=-Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        A_c_right=linspace(2,N-1,N-2);
        A_r_right=linspace(2,N-1,N-2)+1;
        A_0_right=-Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        A_c=[A_c0,A_c_left,A_c_right];
        A_r=[A_r0,A_r_left,A_r_right];
        A_val=[A_0',A_0_left',A_0_right'];
        A=sparse(A_c,A_r,A_val,N,N);
        
        %
        B_c0=linspace(1,N,N);
        B_r0=linspace(1,N,N);
        B_0=ones(N,1);
        B_0(2:N-1)=1-Delta_t/2/d*alpha(T_0(2:N-1))./ct(T_0(2:N-1),Tc,I_0,Ic0)-...
            2*Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        B_c_left=linspace(2,N-1,N-2);
        B_r_left=linspace(2,N-1,N-2)-1;
        B_0_left=Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        B_c_right=linspace(2,N-1,N-2);
        B_r_right=linspace(2,N-1,N-2)+1;
        B_0_right=Delta_t/2/(Delta_x)^2./ct(T_0(2:N-1),Tc,I_0,Ic0).*kappa(T_0(2:N-1),Tc,I_0,Ic0,Rd,d);
        B_c=[B_c0,B_c_left,B_c_right];
        B_r=[B_r0,B_r_left,B_r_right];
        B_val=[B_0',B_0_left',B_0_right'];
        B=sparse(B_c,B_r,B_val,N,N);
        %
        J_0=I_0/w/d;
        C=rho(I_0,Ic0,T_0,Tc,Rd,d)*(J_0)^2*Delta_t./ct(T_0,Tc,I_0,Ic0)+...
            2*Tsub*Delta_t/2/d*alpha(T_0)./ct(T_0,Tc,I_0,Ic0);
        C(1)=0;
        C(N)=0;
        %
        T_ev(:,i)=A\(B*T_0+C);
        
        
        
        R_00=Resistance(T_ini,Tc,I_0,Ic0,Rd,w,Delta_x);
        R_1=Resistance(T_ev(:,i),Tc,I_0,Ic0,Rd,w,Delta_x);
        I_A=ones(2,2);
        I_A(1,1)=2*Lk/Delta_t+Z0+R_1+Rs;
        I_A(1,2)=-1;
        I_A(2,1)=Delta_t/2/Cbt;
        I_B=ones(2,2);
        I_B(1,1)=2*Lk/Delta_t-Z0-R_ev(i-1)-Rs;
        I_B(2,1)=-Delta_t/2/Cbt;
        I_C=[2*(Z0)*Ibias;2*Delta_t/2/Cbt*Ibias];%Z0 with Rs or not
        x0=[I_ev(i-1);Vc_ev(i-1)];
        x1=I_A\(I_B*x0+I_C);
        I_ev(i)=x1(1);
        Vc_ev(i)=x1(2);
        R_ev(i)=R_1;
    end
    T_ev_cleo(:,:,R_num)=T_ev;
    I_ev_cleo(:,:,R_num)=I_ev;
    Vc_ev_cleo(:,:,R_num)=Vc_ev;
    R_ev_cleo(:,:,R_num)=R_ev;
    R_num=R_num+1;
end
%%

%% 

%% CLEO-Temperature
T_label=linspace(1,step_t+2,step_t+2)*Delta_t;
x_label=linspace(-N/2+0.5,N/2-0.5,N+2)*Delta_x;
set(0,'defaultfigurecolor','w')
figure()
t=tiledlayout(1,3,'TileSpacing','compact');

%tiledlayout('flow');
set(gcf,'Position',[100 100 1250 400]);%origin and length width
%set(gca,'Fontsize',20)
%set(gca,'Position',[.13 .17 .60 .74]); 

ax1=nexttile(1,[1,1]);
imagesc(T_label*10^9,x_label*10^6,T_ev_cleo(:,:,1))
% set(gca,'Fontsize',18,'FontName','Times New Roman')
xlim([0,3])
ylim([-0.65,0.65])
text(0.16,-0.5,['{\it R}_s = 0 ',char(937)],'Interpreter','tex','FontSize',18,'Color','w')
set(gca,'Fontsize',18,'FontName','Helvetica')
%ylim([-3,3])
% title('(c)','Units','Normalized','Position',[-0.1,1.1,0])


ax2=nexttile(2,[1,1]);
imagesc(T_label*10^9,x_label*10^6,T_ev_cleo(:,:,2))
% set(gca,'Fontsize',18,'FontName','Times New Roman')
xlim([0,3])
ylim([-0.45,0.45])
text(0.16,-0.5,['{\it R}_s = 150 ',char(937)],'Interpreter','tex','FontSize',18,'Color','w')
set(gca,'Fontsize',18,'FontName','Helvetica')
%ylim([-3,3])
%title('(c)','Units','Normalized','Position',[-0.1,1.1,0])
colormap('jet')


ax3=nexttile(3,[1,1]);
imagesc(T_label*10^9,x_label*10^6,T_ev_cleo(:,:,3))
% set(gca,'Fontsize',18,'FontName','Times New Roman')
xlim([0,3])
text(0.16,-0.5,['{\it R}_{s} = 350 ',char(937)],'Interpreter','Tex','FontSize',18,'Color','w')
ylim([-0.45,0.45])
%title('(c)','Units','Normalized','Position',[-0.1,1.1,0])
linkaxes([ax1,ax2,ax3],'xy')
cb=colorbar('position',[0.925,0.222,0.02,0.702]);
colormap('jet')
xlabel(t,'Time (ns)','Interpreter','tex','FontSize',22)
ylabel(t,['x (',char(956),'m)'],'Interpreter','tex','FontSize',22)%https://unicodelookup.com/#greek
ylabel(cb,'T (K)','Interpreter','tex','FontSize',18)
set(gca,'Fontsize',18,'FontName','Helvetica')
cb.Label.Position(1)=2;
%cb.Layout.Tile='east';


%save images
% filename='Paper-temp-nbtin-20250806-300Rd-7_707K_2.5K';
% saveas(gcf,[filename,'.fig']);
% print(gcf,'-dsvg',[filename,'.svg'])
% exportgraphics(gcf,[filename,'.png'],'Resolution',600);
%% CLEO-Voltage 
T_label=linspace(1,step_t+2,step_t+2)*Delta_t;
x_label=linspace(-N/2+0.5,N/2-0.5,N+2)*Delta_x;
set(0,'defaultfigurecolor','w')
figure()

set(gcf,'Position',[100 100 700 400]);%origin and length width

plot(T_label*10^9,(Ibias-I_ev_cleo(:,:,1))/(Ibias-min(I_ev_cleo(:,:,1))),'DisplayName',['0 ',char(937)],'color',[18/255,80/255,123/255],'LineWidth',2);
tup=T_label(find(I_ev_cleo(:,:,1)==min(I_ev_cleo(:,:,1))));
% tfall_index=find(I_ev_cleo(:,:,1)==Ibias-1/exp(2)*(Ibias-min(I_ev_cleo(:,:,1))));
[~,tfall_index]=min(abs(Ibias-1/exp(1)*(Ibias-min(I_ev_cleo(1:end,:,1)))-I_ev_cleo(100:end,:,1)));
tfall=T_label(tfall_index(end)+99);
t1=tfall-tup;%1.2ns
hold on
plot(T_label*10^9,(Ibias-I_ev_cleo(:,:,2))/(Ibias-min(I_ev_cleo(:,:,2))),'DisplayName',['150 ',char(937)],'color',[250/255,192/255,61/255],'LineWidth',2);
tup2=T_label(find(I_ev_cleo(:,:,2)==min(I_ev_cleo(:,:,2))));
[~,tfall_index2]=min(abs(Ibias-1/exp(1)*(Ibias-min(I_ev_cleo(1:end,:,2)))-I_ev_cleo(100:end,:,2)));
tfall2=T_label(tfall_index2(end)+99);
t2=tfall2-tup2;%0.234ns
hold on
plot(T_label*10^9,(Ibias-I_ev_cleo(:,:,3))/(Ibias-min(I_ev_cleo(:,:,3))),'DisplayName',['350 ',char(937)],'color',[209/255,41/255,32/255],'LineWidth',2);
hold off

xlabel('Time (ns)','Interpreter','tex','FontSize',22)
ylabel('Normalized voltage (a.u.)','Interpreter','tex','FontSize',22)
lgd=legend;
lgd.FontSize=17;
lgd.EdgeColor='k';
lgd.NumColumns=1;
lgd.Location='east';

title(lgd,'{\it R}_s','Interpreter','tex','FontSize',20);
lgd.Box='off';
set(lgd,'Fontsize',18,'Interpreter','tex');
% grid on
box on
ylim([0,1.1])
set(gca,'Fontsize',18,'FontName','Helvetica')

%save images
% filename='Paper-volt-nbtin-20250806-300Rd-7_707K_2.5K';
% saveas(gcf,[filename,'.fig']);
% print(gcf,'-dsvg',[filename,'.svg'])
% exportgraphics(gca,[filename,'.png'],'Resolution',600);
%% fit falling edge
tup=T_label(find(I_ev_cleo(:,:,2)==min(I_ev_cleo(:,:,2))));
% tfall_index=find(I_ev_cleo(:,:,1)==Ibias-1/exp(2)*(Ibias-min(I_ev_cleo(:,:,1))));
[~,tfall_index]=min(abs(Ibias-1/exp(1)*(Ibias-min(I_ev_cleo(1:end,:,2)))-I_ev_cleo(100:end,:,2)));
% [~,tfall_index]=min(abs(Ibias-1/exp(1)*(Ibias-min(I_ev_cleo(1:end,:,1)))-I_ev_cleo(100:end,:,1)));
tfall=T_label(tfall_index(end)+99);

fit_predict_time=T_label(T_label>tup)'*1e9;
fit_predict_eff=(Ibias-I_ev_cleo((T_label>tup),:,2))*50*1000;
% fit_predict_eff=(Ibias-I_ev_cleo((T_label>tup),:,1))*50*1000;
options = optimset('Display', 'off');
initial_guess=[max(fit_predict_eff),0,5];
parameters = lsqcurvefit(@exp_fall_func, initial_guess, fit_predict_time, fit_predict_eff, [], [], options);
fit_fall_time=parameters(3);


%% RSNSPD 
T_label=linspace(1,step_t+2,step_t+2)*Delta_t;
x_label=linspace(-N/2+0.5,N/2-0.5,N+2)*Delta_x;
set(0,'defaultfigurecolor','w')
figure(24)

set(gcf,'Position',[100 100 700 400]);%origin and length width

plot(T_label*10^9,(R_ev_cleo(:,:,1)),'DisplayName',['0 ',char(937)],'color',[18/255,80/255,123/255],'LineWidth',2);
tup=T_label(find(I_ev_cleo(:,:,1)==min(I_ev_cleo(:,:,1))));
[~,tfall_index]=min(abs(Ibias-1/exp(2)*(Ibias-min(I_ev_cleo(1:end,:,1)))-I_ev_cleo(100:end,:,1)));
tfall=T_label(tfall_index(end)+99);
t1=tfall-tup;%1.2ns
hold on
plot(T_label*10^9,(R_ev_cleo(:,:,2)),'DisplayName',['150 ',char(937)],'color',[250/255,192/255,61/255],'LineWidth',2);
tup2=T_label(find(I_ev_cleo(:,:,2)==min(I_ev_cleo(:,:,2))));
[~,tfall_index2]=min(abs(Ibias-1/exp(2)*(Ibias-min(I_ev_cleo(1:end,:,2)))-I_ev_cleo(100:end,:,2)));
tfall2=T_label(tfall_index2(end)+99);
t2=tfall2-tup2;%0.234ns
hold on
plot(T_label*10^9,(R_ev_cleo(:,:,3)),'DisplayName',['300 ',char(937)],'color',[209/255,41/255,32/255],'LineWidth',2);
hold off

xlabel('Time (ns)','Interpreter','tex','FontSize',22)
ylabel(['Resistance (',char(937),')'],'Interpreter','tex','FontSize',22)
lgd=legend;
lgd.FontSize=17;
lgd.EdgeColor='k';
lgd.NumColumns=1;
lgd.Location='northeast';

title(lgd,'{\it R}_s','Interpreter','tex','FontSize',20);
lgd.Box='off';
set(lgd,'Fontsize',18,'Interpreter','tex');
% grid on
box on
% ylim([0,1.1])
xlim([0,3])
set(gca,'Fontsize',18,'FontName','Helvetica')

%save images
% filename='Paper-Rsnspd-nbtin-20250408-300Rd-7_644K';
% saveas(gcf,[filename,'.fig']);
% print(gcf,'-dsvg',[filename,'.svg'])
% exportgraphics(gca,[filename,'.png'],'Resolution',600);

%%
function y = exp_fall_func(x,t)
    A=x(1);%amplitude
%     offset=x(4);%vertical offset
    slope=x(3);%slope of the curve
    shift=x(2);%shift
    y = A.*exp(-(t-shift)./slope);
end
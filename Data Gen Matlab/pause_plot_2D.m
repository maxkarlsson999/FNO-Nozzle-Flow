 
function pause_plot_2D( x, y, data2, str_NameOfNumericalScheme, BC)

    global gamma R;
    % clear previous title

    %fh = findobj( 'Type', 'Figure', 'Name', 'Numerical Solutions' );
    
    delete(findall(gcf,'Tag','SuperTitle'));
   
    strTitle = [str_NameOfNumericalScheme ': t = ' num2str(data2.CurrentTime*1000) 'ms'];

    if contains( BC.type , 'Nozzle')
        global x_1D  data_Anlytical_Quasi1D;
        data1 = data_Anlytical_Quasi1D;
        
        subplot(3,4,1); 
        plot( x_1D, data2.P(:,end/2)/100000 ,  '-ro' ); hold on;
        plot( x_1D, data1.P/100000   ,  'k'  ,'LineWidth' , 2);  hold on;
        xlim([x_1D(1),x_1D(end) ]); ylabel( 'P (Bar)'); title( 'Pressure '); %grid on;
        hold off;
        
        subplot(3,4,2); 
        plot( x_1D, data2.rho(:,end/2),  '-ro' ) ;    hold on; 
        plot( x_1D, data1.rho  ,  'k'  ,'LineWidth' , 2);    hold on; 
        xlim([x_1D(1),x_1D(end) ]); ylabel( 'rho (kg/m^3)'); title( 'Density');% grid on;  
        hold off;

        subplot(3,4,5); 
        plot( x_1D, data2.s(:,end/2) ,  '-ro' ); hold on;
        plot( x_1D, data1.s   ,  'k' ,'LineWidth' , 2 ); hold on;
        xlim([x_1D(1),x_1D(end) ]); ylabel( 's (J/K/kg)'); title( 'Entropy');% grid on;
        hold off;

        subplot(3,4,6); 
        plot( x_1D, data2.T(:,end/2) ,  '-ro' ); hold on;
        plot( x_1D, data1.T   ,  'k' ,'LineWidth' , 2 ); %hold on;
        xlim([x_1D(1),x_1D(end) ]); ylabel( 'T (K)'); title( 'Temperature ');% grid on;
        hold off;
    
        subplot(3,4,9);
        velo_centerline = sqrt(data2.u(:,end/2).^2+data2.v(:,end/2).^2);
        plot( x_1D,         velo_centerline,  '-ro' ); hold on;
        plot( x_1D, data1.u,  'k'  ,'LineWidth' , 2); hold on; 
        xlim([x_1D(1),x_1D(end) ]);xlabel( 'xCen(m)'); ylabel( 'u (m/s)'); title( 'Velocity '); %grid on;
        hold off;

        subplot(3,4,10); 
        %plot( x, data2.M ,  '-ro' ); hold on;
        plot( x_1D, velo_centerline./sqrt( gamma*R*data2.T(:,end/2)  )  ,  '-ro' ); hold on;
        plot( x_1D, data1.M   ,  'k' ,'LineWidth' , 2 ); hold on;
        xlim([x_1D(1),x_1D(end) ]); xlabel( 'xCen(m)'); ylabel( 'M'); title( 'Local Mach number'); %grid on; 
        hold off;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%   2D plot   %%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,4, [3 4] )
        contourf(x,y, (data2.P)/10^5  );     colorbar(); title( ['P [bar] ' ] ) ;   daspect([ 1 1 1])
        subplot(3,4, [7 8])
        contourf(x,y, (data2.T)  );     colorbar();       title( ['T [k]  ' ] ) ;  daspect([ 1 1 1])
        subplot(3,4, [11 12 ])
        contourf(x,y, sqrt(data2.u.^2 + data2.v.^2)   );     colorbar();    title( ['|velo| [m/s]  '] ) ;    daspect([ 1 1 1])
    elseif contains(BC.type, 'Dumbbell')
        subplot(3,1, 1 )
        contourf(x,y, (data2.P)/10^5  );     colorbar(); title( ['P [bar] ' ] ) ;   daspect([ 1 1 1])
        subplot(3,1, 2)
        contourf(x,y, (data2.T)  );     colorbar();       title( ['T [k]  ' ] ) ;  daspect([ 1 1 1])
        subplot(3,1, 3)
        contourf(x,y, sqrt(data2.u.^2 + data2.v.^2)   );     colorbar();    title( ['|velo| [m/s]  '] ) ;    daspect([ 1 1 1])
        % colorbar('southoutside') ; %hold on     ;          %quiver(x,y,data.u ,data.v);    
    else
        subplot(2,2, 1 )
        contourf(x,y, (data2.P)/10^5  );     colorbar(); title( ['P [bar] ' ] ) ;   daspect([ 1 1 1])
        subplot(2,2, 2)
        contourf(x,y, (data2.T)  );     colorbar();       title( ['T [k]  ' ] ) ;  daspect([ 1 1 1])
        subplot(2,2, 3)
        contourf(x,y, data2.rho   );     colorbar();    title( ['rho [kg/m^3]  '] ) ;    daspect([ 1 1 1])
        subplot(2,2, 4)
        contourf(x,y, sqrt(data2.u.^2 + data2.v.^2)   );     colorbar();    title( ['|velo| [m/s]  '] ) ;    daspect([ 1 1 1])
        % colorbar('southoutside') ; %hold on     ;          %quiver(x,y,data.u ,data.v);    
        %%%%%%%%%%%%%%%%%%%%%%%%%%   2D plot   %%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    annotation( 'textbox', [0 0.9 1 0.1],  'String', strTitle, 'EdgeColor', 'none',  'FontSize',12,  'HorizontalAlignment', 'center','Tag' , 'SuperTitle');
    

    pause (0.01);
end


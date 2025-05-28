function plotNozzle(x,A)
   f=figure(3); 
   %set(f, 'MenuBar', 'none');    set(f, 'ToolBar', 'none'); 
   set(f, 'Position', [0, 800, 500, 250]);  
   plot(x, A, 'k'); 
   hold on; 
   plot(x, -A,'k');
   title('Nozzle shape');  
end

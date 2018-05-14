maxi = max(sourceInt.pow);
indexymax = find(sourceInt.pow > 0.80*maxi);
pozicemax = sourceInt.pos(indexymax,:);

mini = min(sourceInt.pow);
indexymin = find(sourceInt.pow < 1.02 *mini);
pozicemin = sourceInt.pos(indexymin,:);
zdrojpozice = sourcemodel.pos(zdroj,:);


figure;
hold on;
scatter3(pozicemax(:,1),pozicemax(:,2),pozicemax(:,3),...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor','b');
scatter3(zdrojpozice(1),zdrojpozice(2),zdrojpozice(3),...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor','r');
   scatter3(pozicemin(:,1),pozicemin(:,2),pozicemin(:,3),...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor','g');
 
   scatter3(gridinside(:,1),gridinside(:,2),gridinside(:,3),...
       'MarkerEdgeColor','k');

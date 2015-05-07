hFig = figure('Menubar','none');
table = uitable('Parent',hFig);
data = {'sdf', [12; 23], false};
% data = {[], [], [];[],[],[]};
set(table, 'Units','normalized', 'Data', data,...
                'ColumnName', {'char', 'numeric', 'logical'}, 'ColumnFormat', {'char', 'numeric', 'logical'},...
                'ColumnEditable', [true, true, true]);
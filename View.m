classdef View
    
    properties
        handles;
        tableRows;
    end
    
    methods
        
        function this = View(M)
            %VIEW  a GUI representation of the signal model
            handles_ = View.initGUI();
            this.tableRows = this.formIsotermTable(M, handles_.tabOut);
            
            plot(M.Cr, M.Ar,'greeno', 'LineWidth', 3, 'Parent', handles_.axes);
            x = linspace(0, max(M.Cr), 1000);
            hold on;
            for i = 1:length(M.isoterms)
                isoterm = M.isoterms{i};
                if(~isempty(isoterm) && isprop(isoterm, 'isotermResult'))
                    y = isoterm.isotermResult(x);
                    plot(x, y, 'Parent', handles_.axes);
                    display(isoterm.isotermResult);
                end
            end
            legend('off');
            
            this.handles = handles_;
        end
        
        
        
    end
    
    methods(Static)
        function tableRows = formIsotermTable(M, isotermTable)
            N = length(M.isoterms);
            tableRows = cell(N, 1);
            tableRowsData = cell(N, IsotermTableRow.size);
            for i = 1:N
                if isempty(M.isoterms{i})
                     tableRows{i} = IsotermTableRow(M.isotermTypes{i});
                else
                     tableRows{i} = IsotermTableRow(M.isoterms{i});
                end
                
                tableRowsData(i, :) = tableRows{i}.data;
            end
            set(isotermTable, 'Units','normalized', 'Data', tableRowsData,...
                'ColumnName', IsotermTableRow.columnName, 'ColumnFormat', IsotermTableRow.columnFormat,...
                'ColumnEditable', IsotermTableRow.columnEditable);
            columnWidth = cell(1, N);
            columnWidth(:) = {'auto'};
            set(isotermTable,'ColumnWidth',columnWidth);
        end
        
        function handles = initGUI()
            % load FIG file (its really a MAT-file)
            hFig = hgload('View.fig');
            % extract handles to GUI components
            hAx = findobj(hFig, 'tag','axes');
            hCheckConfInt = findobj(hFig, 'tag','checkConfInt');
            hEditConfInt = findobj(hFig, 'tag','editConfInt');
            hButRecalc = findobj(hFig, 'tag','butRecalc');
            hTabIn = findobj(hFig, 'tag','tabIn');
            hTabOut = findobj(hFig, 'tag','tabOut');
            
            hButLoadData = findobj(hFig, 'tag','butLoadData');
            
            % return a structure of GUI handles
            handles = struct('fig',hFig, 'axes',hAx, 'checkConfInt',hCheckConfInt, ...
                'editConfInt',hEditConfInt, 'butLoadData',hButLoadData, 'butRecalc',...
                hButRecalc, 'tabIn', hTabIn, 'tabOut', hTabOut);
        end
    end
end

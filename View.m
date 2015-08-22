classdef View < handle
    
    properties
        hGUI;
        hModel;
    end
    
    properties(SetObservable = true)
        tableRows;%?need
        axisHandles = []
        hDataAxes
    end
    
    methods
        
        function this = View(Model)
            %VIEW  a GUI representation of the signal model
            handles = View.initGUI();

            this.tableRows = this.formIsotermTable(Model, handles.tabOut);
            isotermsListener = event.proplistener(Model, findprop(Model,'isoterms'), 'PostSet',...
                @(o, e)this.update(o, e));
            setappdata(handles.fig, 'proplistener',isotermsListener);
            
            modelDataListener = event.proplistener(Model, findprop(Model,'data'), 'PostSet',...
                @(o, e)this.updateData());
            setappdata(handles.fig, 'proplistener',modelDataListener);
            
            this.hGUI = handles;
            this.hModel = Model;
        end
        
        function update(this, o, e)
                isotermId = this.hModel.lastIsoInd;
                if(isotermId == -1)
                    this.formIsotermTable(this.hModel, this.hGUI.tabOut);
                    return;
                end
                isoterm = this.hModel.isoterms{isotermId};
                if(isempty(isoterm) || isempty(isoterm.isotermResult))
                    isoterm = this.hModel.isotermTypes{isotermId};
    %                 delete(this.axisHandles(isotermId));
                end   
                table = this.hGUI.tabOut;
                %table
                tableRow = IsotermTableRow(isoterm);
                this.tableRows{isotermId} = tableRow;

                tableRowsData = get(table, 'Data');
                tableRowsData(isotermId, :) = tableRow.data;
                set(table, 'Data', tableRowsData);
            
                %ploting
                if(~isempty(isoterm) && isprop(isoterm, 'isotermResult'))
                    h = plot(isoterm.isotermResult);
                    this.axisHandles(isotermId) = h;
                end
                legend('off');
        end
        
        function updateData(this)
            data = this.hModel.data;
            delete(this.hDataAxes);
            cla;
            handles = this.hGUI;
            axes(handles.axes);
            data(end, :) = [0, 0];
            Cr = data(:, 1);
            Ar = data(:, 2);
            this.hDataAxes = plot(Cr, Ar, 'greeno', 'LineWidth', 3);
            this.axisHandles = zeros(length(Cr), 1);
            hold on;
        end
        
        function isotermIdList = getCheckedRows(this)
            tableOut = this.hGUI.tabOut;
            dataOut = get(tableOut, 'Data');
            isShownColumn = dataOut(:, IsotermTableRow.columnShow);
            isotermIdList = find(cell2mat(isShownColumn) == 1)';
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
            set(isotermTable, 'Data', tableRowsData,...
                'ColumnName', IsotermTableRow.columnName, 'ColumnFormat', IsotermTableRow.columnFormat,...
                'ColumnEditable', IsotermTableRow.columnEditable);
            
            columnWidth = max(cellfun('length', tableRowsData));
            columnWidth(columnWidth < 6) = 0;
            columnWidth = columnWidth * 7;
            columnWidthCell = num2cell(columnWidth);
            for i = 1:length(columnWidthCell)
                if(columnWidthCell{i} == 0)
                    columnWidthCell{i} = 'auto';
                end
            end
            set(isotermTable,'ColumnWidth',columnWidthCell);
        end
        
        
        function handles = initGUI()
            % load FIG file (its really a MAT-file)
            hFig = hgload('View.fig');
            handles = struct('fig',hFig);
            tags = {'axes', 'checkConfInt', 'editConfInt', 'editConfInt', 'butRecalc', 'tabIn',...
                'tabOut', 'butPlusRow', 'butMinusRow', 'butLoadData', 'menuTabIn', 'contextPaste'};
            % extract handles to GUI components
            for tag = tags
                 tag = char(tag);
                 hTag = findobj(hFig, 'tag', tag);
                 handles.(tag) = hTag;
            end
            
            set(handles.tabIn, 'uicontextmenu', handles.menuTabIn);
        end
        
    end
end

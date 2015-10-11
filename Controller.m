classdef Controller < handle
    properties
        model;
        view;
        data = [];
    end
    methods
        function this = Controller()
            model = Model();
            view = View(model);
            set(view.hGUI.contextPaste, 'Callback',{@Controller.copyExcel,view.hGUI.tabIn});
            set(view.hGUI.butPlusRow, 'Callback',{@Controller.onClickPlusRow,view.hGUI.tabIn});
            set(view.hGUI.butMinusRow, 'Callback',{@Controller.onClickMinusRow,view.hGUI.tabIn});
            set(view.hGUI.tabOut, 'CellEditCallback', @this.onTableEdit);
            set(view.hGUI.butRecalc, 'Callback', @this.onClickRecalculate);
            set(view.hGUI.cbP1, 'Callback', @this.onCheckP1);
            set(view.hGUI.cbConfInt, 'Callback', @this.onCheckConfInt);
            set(view.hGUI.edConfInt, 'Callback', @this.onFinishConfInt);
            %2,4DNP
            Cr = [
                0.5
                1
                2
                4
                8
                18
                46
                50
                80
                90
                ];
            
            Ar = [
                83.25
                99.8
                124.5
                165.3333333
                246
                337.4
                385.9
                450
                504
                546.6666667
                ];
            data = [Cr Ar];
            set(view.hGUI.tabIn, 'data', data);
            this.model = model;
            this.view = view;
        end
        function onFinishConfInt(this, edConfInt, event)
            tableRows = this.view.tableRows;
%             isotermsId = this.view.getCheckedRows();
            
            isotermIdList = this.view.getAllCalculatedIsoterms();
            tabOutDataRows = get(this.view.hGUI.tabOut, 'Data');
            confInts = tabOutDataRows(:, IsotermTableRow.columnConfInt);
            %this.view.hGUI.edConfInt
            confidenceLevel = str2double(get(edConfInt,'String')) / 100;
            for id = isotermIdList
                confInt = confint(this.model.isoterms{id}.isotermResult, confidenceLevel)';
                confInt(confInt<0) = 0;
                confInts{id} =  mat2str(round(confInt * 100) / 100);
            end
            tabOutDataRows(:, IsotermTableRow.columnConfInt) = confInts;
            set(this.view.hGUI.tabOut, 'Data', tabOutDataRows)
            
            isCheckedConfInt =  get(this.view.hGUI.cbConfInt, 'Value');
            for id = isotermIdList
                hConfInt = tableRows{id}.hConfIntAxes;
                if(ishghandle(hConfInt))
                    delete(hConfInt);
                end
                this.view.tableRows{id}.hConfIntAxes = this.view.plotConfInt(this.model.isoterms{id});
                this.view.tableRows{id}.updatePlot(isCheckedConfInt);
            end
        end
        
        function onCheckConfInt(this, cbConfInt, event)
            isChecked =  get(this.view.hGUI.cbConfInt, 'Value');
            
            chekedIsotermsID = this.view.getCheckedRows();
            for isotermId = chekedIsotermsID
                if(isChecked)
                    isoterm = this.model.isoterms{isotermId};
                    if ~isempty(isoterm)
                        hConfPlot = this.view.plotConfInt(isoterm);
                        this.view.tableRows{isotermId}.hConfIntAxes = hConfPlot;
                    end
                else
                    hConfLine =  this.view.tableRows{isotermId}.hConfIntAxes;
                    if(ishghandle(hConfLine)) 
                        delete(hConfLine);
                    end;
                end
            end
        end
        
        function onCheckP1(this, cbP1, event)
            isChecked =  get(cbP1, 'Value');
            hEdP1 = this.view.hGUI.edP1;
            if(isChecked)
               set(hEdP1, 'Visible', 'off');
            else
               set(hEdP1, 'Visible', 'on');
            end
        end
        
        function onClickRecalculate(this, butRecalc, event)
            isChecked =  get(this.view.hGUI.cbP1, 'Value');
            if(isChecked)
                this.model.parameters = [];
            else
                this.model.parameters = str2num(get(this.view.hGUI.edP1, 'String'));
            end
            isotermsIdList = this.view.getCheckedRows();
            tableIn = this.view.hGUI.tabIn;
            dataIn = get(tableIn, 'Data');
            if(~isequal(this.data, dataIn))
                this.model.data = dataIn;
                this.model.clear();
            end
            this.data = dataIn;
            
            table = this.view.hGUI.tabOut;
            tableRowsData = get(table, 'Data');
            for i = isotermsIdList
                tableRowsData{i, IsotermTableRow.columnShow} = true;
            end
            set(table, 'Data', tableRowsData);
            
            this.model.calculate(isotermsIdList);
            View.fitTableWidth(table);
        end

        function onTableEdit(this, table, eventData)
            row = eventData.Indices(1);
            col = eventData.Indices(2);
            dataOut = get(table, 'Data');
           
            if(col == IsotermTableRow.columnShow)
                tableRow = this.view.tableRows{row};
                tableRow.data{IsotermTableRow.columnShow} = dataOut{row, col}; 
                isCheckedConfInt =  get(this.view.hGUI.cbConfInt, 'Value');
                tableRow.updatePlot(isCheckedConfInt);
                this.view.tableRows{row} = tableRow;
            end
        end
    end
    
    methods(Static = true)
        function copyExcel(o, e, table)
            str = clipboard('paste');
            str = strrep(str,',','.');
            import = str2num(str);
            set(table, 'Data', import);
        end
        
        function onClickPlusRow(o, e, table)
            data = get(table, 'Data');
            data(end+1,:) = [0, 0];
            set(table, 'Data', data);
        end
        
        function onClickMinusRow(o, e, table)
            data = get(table, 'Data');
            set(table, 'Data', data(1:end-1,:));
        end
        
    end
end

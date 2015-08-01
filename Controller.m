classdef Controller < handle
    properties
        model;
        view;
    end
    methods
        function this = Controller()
            model = Model();
            view = View(model);
            set(view.hGUI.contextPaste, 'Callback',{@Controller.copyExcel,view.hGUI.tabIn})
            set(view.hGUI.butPlusRow, 'Callback',{@Controller.onClickPlusRow,view.hGUI.tabIn})
            set(view.hGUI.butMinusRow, 'Callback',{@Controller.onClickMinusRow,view.hGUI.tabIn})
            set(view.hGUI.butLoadData, 'Callback',@this.onClickLoadData)
           
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
        
        function onClickLoadData(this, o, e)
            tableIn = this.view.hGUI.tabIn;
            tableOut = this.view.hGUI.tabOut;
            dataIn = get(tableIn, 'Data');
            dataOut = get(tableOut, 'Data');
            isShown = dataOut(:, IsotermTableRow.columnShow);
            this.model.data = dataIn;
            
            ar = find(cell2mat(isShown) == 1);
            this.model.calculate(ar');
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
            data(end+1,:) = {[]};
            set(table, 'Data', data);
        end
        
        function onClickMinusRow(o, e, table)
            data = get(table, 'Data');
            set(table, 'Data', data(1:end-1,:));
        end
        
    end
end

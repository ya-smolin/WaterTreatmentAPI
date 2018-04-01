classdef DFun < handle
    properties
        grid % HxMxN
        dim
        gridSize % H, M, N
        varSpace % 1:H 1:M 2:N
    end
    
    methods
        %(x, r, t) or (x, u, t)
        function o = DFun(gridSize)
            o.dim = length(gridSize);
            o.grid = zeros(gridSize);
            o.gridSize = gridSize;
            o.varSpace = ones(o.dim, 2);
            o.varSpace(:,2) = o.gridSize;
        end
        
        function condInitial(o, funX0)
            o.condBoundLow(funX0, o.dim)
        end
        
        function condBoundLow(o, fun0T, num)
            o.boundaryCondConst(fun0T, 1, num)
            o.varSpace(num, 1) = 2;
        end
        
        function condBoundUp(o, funHT, num)
            o.boundaryCondConst(funHT, o.gridSize(num), num)
            o.varSpace(num, 2) = o.gridSize(num) - 1;
        end
        
        function boundaryCondConst(o, value, placeInd, num)
            if o.dim == 2
                if num == 1
                    o.grid(placeInd, :) = value;
                elseif num == 2
                    o.grid(:, placeInd) = value;
                end
            elseif o.dim == 3
                if num == 1
                    o.grid(placeInd, :, :) = value;
                elseif num == 2
                    o.grid(:, placeInd, :) = value;
                elseif num == 3
                    o.grid(:, :, placeInd) = value;
                else
                    error("num have to be less or equal dim");
                end
            else
                error("current support only 2D and 3D. implement here if necessary");
            end
        end
        
        function sizeGridVar = sizeGridVar(o)
            sizeGridVar = prod(o.varSize(1:o.dim-1));
        end
        
        function varSize = varSize(o, indexes)
            varSize = diff(o.varSpace') + 1;
            if nargin > 1
                varSize = varSize(indexes);
            end
        end
        
        function varDia = varDia(o, num)
            varDia = o.varSpace(num, 1):o.varSpace(num, 2);
        end
        
        function unpack(o, X, k)
            if o.dim == 2
                o.grid(o.varDia(1), k) = X;
            elseif o.dim == 3
                o.grid(o.varDia(1), o.varDia(2), k) = reshape(X, o.varSize(1), o.varSize(2), 1);
            else
                error("implemented only for 2D and 3D. implement here if necessary");
            end
        end
        
        function pass = unitTest3D(o)
            testObj = DFun([3, 4, 5]);
            MOCK_X = 22;
            MOCK_Y = 33;
            MOCK_T = 11;
            testObj.condBoundLow(MOCK_X, 1);
            testObj.condBoundUp(MOCK_Y, 2);
            testObj.condInitial(MOCK_T);
            
            [ansStr, pass] = equalEps(testObj.varSpace, [2 3; 1 3; 2 5]);
            disp(['varSpace: ' ansStr]);
            [ansStr, pass1] =  equalEps(testObj.varDia(3), 2:5);
            disp(['varDia: ' ansStr]);
            [ansStr, pass2] =  equalEps(testObj.varSize(1:3), [2, 3, 4]);
            disp(['varSize: ' ansStr]);
            [ansStr, pass3] =  equalEps(testObj.sizeGridVar(), 6);
            disp(['sizeGridVar: ' ansStr]);
            
            gridTest = zeros(3, 4);
            gridTest(1,:) = MOCK_X;
            gridTest(:,4) = MOCK_Y;
            gridReal = squeeze(testObj.grid(:, :, 2));
            [ansStr, pass4] =  equalEps(gridTest, gridReal);
            disp(['grid: ' ansStr]);
            pass = pass && pass1 && pass2 && pass3 && pass4;
        end
        
        function pass = unitTest2D(o)
            testObj = DFun([3, 4]);
            MOCK_X = 22;
            MOCK_T = 11;
            testObj.condBoundLow(MOCK_X, 1);
            testObj.condInitial(MOCK_T);
            
            [ansStr, pass] = equalEps(testObj.varSpace, [2 3; 2 4]);
            disp(['varSpace: ' ansStr]);
            [ansStr, pass1] =  equalEps(testObj.varDia(2), 2:4);
            disp(['varDia: ' ansStr]);
            [ansStr, pass2] =  equalEps(testObj.varSize(), [2, 3]);
            disp(['varSize: ' ansStr]);
            [ansStr, pass3] =  equalEps(testObj.sizeGridVar(), 2);
            disp(['sizeGridVar: ' ansStr]);
            
            gridTest = zeros(3, 1);
            gridTest(1) = MOCK_X;
            gridReal = testObj.grid(:, 2);
            [ansStr, pass4] =  equalEps(gridTest, gridReal);
            disp(['grid: ' ansStr]);
            pass = pass && pass1 && pass2 && pass3 && pass4;
        end
        
    end
end


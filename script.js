const TOLERANCE = 1e-9; 

// ====== Fraction Converter ======
function formatNumber(num) {
    if (typeof num === 'string') return num; // Return symbol if string
    if (Math.abs(num) < TOLERANCE) return "0";
    if (Math.abs(num - Math.round(num)) < TOLERANCE) {
        return Math.round(num).toString();
    }
    for (let d = 2; d <= 100; d++) {
        let n = num * d;
        if (Math.abs(n - Math.round(n)) < 1e-5) {
            return `${Math.round(n)}/${d}`;
        }
    }
    return Number(num.toFixed(2)).toString();
}

function getOperationText(row1, row2, factor, opType) {
    if (Math.abs(factor) < TOLERANCE && opType !== 'swap') return ''; 
    
    if (opType === 'swap') {
        return `Swap R${row1} <-> R${row2}`;
    }

    if (opType === 'multiply') {
        let factorStr = formatNumber(factor);
        if (factor < 0 && factorStr.includes('/')) factorStr = `(${factorStr})`;
        return `R${row1} = ${factorStr} * R${row1} (Make Pivot 1)`;
    }

    let factorAbs = Math.abs(factor);
    let sign = factor > 0 ? '-' : '+';
    let factorString = formatNumber(factorAbs);
    if (Math.abs(factorAbs - 1) < TOLERANCE) factorString = ''; 

    return `R${row1} = R${row1} ${sign} ${factorString} R${row2}`;
}

// Helper to render a standalone matrix grid
function renderMatrixGrid(matrixData, isCalculation = false) {
    const grid = document.createElement('div');
    grid.className = 'matrix-display';
    let cols = matrixData[0].length;
    grid.style.gridTemplateColumns = `repeat(${cols}, auto)`; 

    for (let i = 0; i < matrixData.length; i++) {
        for (let j = 0; j < cols; j++) {
            const cell = document.createElement('div');
            cell.className = 'matrix-cell';
            if (isCalculation) cell.classList.add('calc-cell');
            
            let val = matrixData[i][j];
            // If it's a number, format it. If string (like "x1"), keep it.
            cell.textContent = typeof val === 'number' ? formatNumber(val) : val;
            grid.appendChild(cell);
        }
    }
    return grid;
}

// Updated displayMatrix to handle custom separator for [A | I]
function displayMatrix(matrix, M, N, container, separatorCol = -1) {
    const grid = document.createElement('div');
    grid.className = 'matrix-display';
    grid.style.gridTemplateColumns = `repeat(${matrix[0].length}, auto)`; 

    for (let i = 0; i < M; i++) {
        for (let j = 0; j < matrix[0].length; j++) {
            const cell = document.createElement('div');
            cell.textContent = formatNumber(matrix[i][j]); 
            cell.className = 'matrix-cell';

            // Logic to draw the vertical bar
            if (separatorCol !== -1) {
                // For Inverse: separator is after A (at index M-1, so border right)
                if (j === separatorCol - 1) cell.style.borderRight = '3px solid #000';
            } else if (j === N) {
                // Default for Gaussian: separator is before b (last column)
                cell.style.borderLeft = '3px solid #000';
            }
            
            grid.appendChild(cell);
        }
    }
    container.appendChild(grid);
}

function displayStep(operation, matrix, M, N, separatorCol = -1) {
    if (operation === "") return;
    const stepsContainer = document.getElementById('steps-container');
    const stepDiv = document.createElement('div');
    stepDiv.className = 'step-item';
    
    const operationText = document.createElement('span');
    operationText.className = 'operation-text';
    operationText.textContent = operation;
    stepDiv.appendChild(operationText);
    
    displayMatrix(matrix, M, N, stepDiv, separatorCol);
    stepsContainer.appendChild(stepDiv);
}

function buildMatrix() {
    const M = parseInt(document.getElementById('num-rows').value); 
    const TotalCols = parseInt(document.getElementById('total-cols').value);
    const N = TotalCols - 1; 
    
    const container = document.getElementById('matrix-container');
    container.innerHTML = ''; 

    if (M < 1 || N < 1 || M > 10 || TotalCols > 11) {
        alert("Please ensure rows (1-10) and total columns (2-11) are within limits.");
        return;
    }

    for (let i = 0; i < M; i++) {
        const rowDiv = document.createElement('div');
        rowDiv.className = 'matrix-row';

        for (let j = 0; j < N; j++) {
            const input = document.createElement('input');
            input.type = 'number';
            input.id = `A-${i}-${j}`;
            input.placeholder = `A[${i+1}, ${j+1}]`;
            rowDiv.appendChild(input);
        }

        const bInput = document.createElement('input');
        bInput.type = 'number';
        bInput.id = `b-${i}`;
        bInput.placeholder = `b[${i+1}]`;
        bInput.style.borderLeft = '3px solid #000'; 
        rowDiv.appendChild(bInput);
        container.appendChild(rowDiv);
    }
    document.getElementById('steps-container').innerHTML = '';
}

function getMatrix(M, N) {
    const matrix = [];
    let allFieldsFilled = true;
    for (let i = 0; i < M; i++) {
        const row = [];
        for (let j = 0; j < N; j++) {
            const inputElement = document.getElementById(`A-${i}-${j}`);
            const value = parseFloat(inputElement ? inputElement.value : NaN);
            if (isNaN(value)) allFieldsFilled = false;
            row.push(value);
        }
        const bElement = document.getElementById(`b-${i}`);
        const bValue = parseFloat(bElement ? bElement.value : NaN);
        if (isNaN(bValue)) allFieldsFilled = false;
        row.push(bValue);
        matrix.push(row);
    }
    if (!allFieldsFilled) return null;
    return matrix;
}

// ====== Core Elimination Logic ======
function forwardElimination(M, colsCount, matrix, stepsContainer, isInverseMode = false) {
    let currentPivotRow = 0;
    let separatorCol = isInverseMode ? M : (colsCount - 1); 
    let pivotLimit = isInverseMode ? M : Math.min(M, colsCount);

    for (let j = 0; j < pivotLimit && currentPivotRow < M; j++) { 

        if (Math.abs(matrix[currentPivotRow][j]) < TOLERANCE) {
            let swapRow = -1;
            for (let i = currentPivotRow + 1; i < M; i++) {
                if (Math.abs(matrix[i][j]) > TOLERANCE) {
                    swapRow = i;
                    break;
                }
            }
            if (swapRow !== -1) {
                [matrix[currentPivotRow], matrix[swapRow]] = [matrix[swapRow], matrix[currentPivotRow]];
                displayStep(getOperationText(currentPivotRow + 1, swapRow + 1, 0, 'swap'), matrix, M, colsCount - 1, separatorCol);
            } else {
                // Found a zero column (pivot not found)
                if (isInverseMode) return "Reason: Column " + (j+1) + " has no pivot (all zeros below diagonal)."; 
                continue; 
            }
        }

        let pivotVal = matrix[currentPivotRow][j];
        // Check if pivot is zero after potential swap (Linear Dependence)
        if (isInverseMode && Math.abs(pivotVal) < TOLERANCE) return "Reason: Row " + (currentPivotRow+1) + " is a linear combination of other rows (Zero Pivot).";

        if (Math.abs(pivotVal - 1) > TOLERANCE) {
            let factor = 1 / pivotVal;
            let operation = getOperationText(currentPivotRow + 1, 0, factor, 'multiply');
            for (let k = j; k < colsCount; k++) {
                matrix[currentPivotRow][k] = matrix[currentPivotRow][k] * factor;
            }
            displayStep(operation, matrix, M, colsCount - 1, separatorCol);
        }

        for (let i = currentPivotRow + 1; i < M; i++) {
            let factor = matrix[i][j]; 
            if (Math.abs(factor) < TOLERANCE) continue; 
            let operation = getOperationText(i + 1, currentPivotRow + 1, factor, 'subtract');
            for (let k = j; k < colsCount; k++) {
                matrix[i][k] = matrix[i][k] - factor * matrix[currentPivotRow][k];
            }
            displayStep(operation, matrix, M, colsCount - 1, separatorCol);
        }
        currentPivotRow++;
    }
    return currentPivotRow; 
}

function backwardElimination(M, colsCount, matrix, stepsContainer, isInverseMode = false) {
        let separatorCol = isInverseMode ? M : (colsCount - 1);
        let pivotLimit = isInverseMode ? M : colsCount - 1;

    let pivotCols = [];
    for (let i = 0; i < M; i++) {
        for (let j = 0; j < pivotLimit; j++) {
            if (Math.abs(matrix[i][j] - 1) < TOLERANCE) {
                pivotCols.push({row: i, col: j});
                break;
            }
        }
    }

    for (const {row: i, col: j} of pivotCols) {
        for (let k = i - 1; k >= 0; k--) {
            let factor = matrix[k][j]; 
            if (Math.abs(factor) < TOLERANCE) continue;
            let operation = getOperationText(k + 1, i + 1, factor, 'subtract');
            for (let l = j; l < colsCount; l++) {
                matrix[k][l] = matrix[k][l] - factor * matrix[i][l];
            }
            displayStep(operation, matrix, M, colsCount - 1, separatorCol);
        }
    }
}

// ====== Enhanced Inverse Matrix Method ======
function solveInverseMethod(M, N, originalMatrix, stepsContainer) {
    // 1. Separate A and b
    let A = [];
    let bVec = [];
    let xVec = [];

    for(let i=0; i<M; i++){
        A.push(originalMatrix[i].slice(0, N));
        bVec.push([originalMatrix[i][N]]); // as 2D column vector
        xVec.push([`x${i+1}`]); // as 2D column vector
    }

    // === 1. Display System Representation (A * X = B) ===
    stepsContainer.innerHTML += '<h3>1. System Representation (Matrix Form):</h3>';
    const sysContainer = document.createElement('div');
    sysContainer.className = 'equation-container';

    // Wrap each matrix with a label
    const wrapMatrix = (grid, label) => {
        const wrapper = document.createElement('div');
        wrapper.className = 'matrix-wrapper';
        const labelDiv = document.createElement('div');
        labelDiv.className = 'matrix-label';
        labelDiv.textContent = label;
        wrapper.appendChild(labelDiv);
        wrapper.appendChild(grid);
        return wrapper;
    };

    const aGrid = renderMatrixGrid(A);
    const xGrid = renderMatrixGrid(xVec);
    const bGrid = renderMatrixGrid(bVec);

    sysContainer.appendChild(wrapMatrix(aGrid, 'Matrix A (Coefficient)'));
    
    const dotOp = document.createElement('div');
    dotOp.className = 'operator';
    dotOp.innerHTML = '&times;';
    sysContainer.appendChild(dotOp);

    sysContainer.appendChild(wrapMatrix(xGrid, 'Vector X'));

    const eqOp = document.createElement('div');
    eqOp.className = 'operator';
    eqOp.innerHTML = '=';
    sysContainer.appendChild(eqOp);

    sysContainer.appendChild(wrapMatrix(bGrid, 'Vector B'));
    
    stepsContainer.appendChild(sysContainer);

    // 2. Build Augmented Matrix [A | I]
    let augMatrix = [];
    for(let i=0; i<M; i++){
        let row = [...A[i]]; 
        for(let k=0; k<M; k++){
            row.push(i === k ? 1 : 0);
        }
        augMatrix.push(row);
    }

    stepsContainer.innerHTML += '<h3>2. Form Augmented Matrix [ A (Coefficient) | I ]:</h3>';
    displayMatrix(augMatrix, M, N, stepsContainer, M); 

    // 3. Apply Gauss-Jordan
    stepsContainer.innerHTML += '<h3>3. Apply Row Operations (Gauss-Jordan) to get Identity on Left:</h3>';
    let pivotResult = forwardElimination(M, 2*M, augMatrix, stepsContainer, true);
    
    // Handle Specific Singular Reasons
    if (typeof pivotResult === 'string') {
        stepsContainer.innerHTML += `<p class="final-result" style="color: red;">Matrix is Singular. ${pivotResult} Cannot calculate Inverse.</p>`;
        return;
    }
    if (pivotResult === -1) {
            // Fallback for generic -1 if returned
            stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Matrix is Singular. Cannot calculate Inverse.</p>';
            return;
    }

    backwardElimination(M, 2*M, augMatrix, stepsContainer, true);

    // Verify Identity
    let isIdentity = true;
    for(let i=0; i<M; i++){
        if(Math.abs(augMatrix[i][i] - 1) > TOLERANCE) isIdentity = false;
    }

    if (!isIdentity) {
            stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Matrix is Singular. Reduced form did not yield Identity matrix.</p>';
            return;
    }

    // Extract Inverse
    let inverseMatrix = [];
    for(let i=0; i<M; i++){
        inverseMatrix.push(augMatrix[i].slice(M, 2*M));
    }

    stepsContainer.innerHTML += '<h3>4. Inverse Matrix (A⁻¹) Found:</h3>';
    displayMatrix(inverseMatrix, M, M, stepsContainer);

    // === 4. Multiplication Steps Visualization ===
    stepsContainer.innerHTML += '<h3>5. Calculate X = A⁻¹ * B:</h3>';
    
    // Container for multiplication equation
    const multContainer = document.createElement('div');
    multContainer.className = 'equation-container';
    
    const invGrid = renderMatrixGrid(inverseMatrix);
    const bGridCalc = renderMatrixGrid(bVec);

    multContainer.appendChild(wrapMatrix(invGrid, 'A⁻¹'));

    const multOp = document.createElement('div');
    multOp.className = 'operator';
    multOp.innerHTML = '&times;';
    multContainer.appendChild(multOp);

    multContainer.appendChild(wrapMatrix(bGridCalc, 'B'));
    
    stepsContainer.appendChild(multContainer);

    // Calculation Matrix Steps
    stepsContainer.innerHTML += '<p class="operation-text">Multiplication Details:</p>';
    let calcMatrixData = [];
    let X = new Array(M).fill(0);

    for(let i=0; i<M; i++){
        let rowSum = 0;
        let stepStr = "";
        let parts = [];
        for(let j=0; j<M; j++){
            let val = inverseMatrix[i][j];
            let bVal = bVec[j][0];
            rowSum += val * bVal;
            
            if (Math.abs(val) > TOLERANCE) {
                parts.push(`(${formatNumber(val)} × ${formatNumber(bVal)})`);
            }
        }
        if (parts.length === 0) parts.push("0");
        stepStr = parts.join(" + ");
        calcMatrixData.push([stepStr]); // Column vector of strings
        X[i] = rowSum;
    }

    // Display calculation grid
    const calcGrid = renderMatrixGrid(calcMatrixData, true);
    stepsContainer.appendChild(calcGrid);

    // Arrow down
    const arrowDiv = document.createElement('div');
    arrowDiv.style.textAlign = 'left';
    arrowDiv.style.padding = '10px 20px';
    arrowDiv.style.fontSize = '1.5em';
    arrowDiv.innerHTML = '&#8595;'; // Down arrow
    stepsContainer.appendChild(arrowDiv);

    // Final Result Matrix
    const resultMatrixData = X.map(val => [val]);
    const resGrid = renderMatrixGrid(resultMatrixData);
    stepsContainer.appendChild(resGrid);


    // 7. Final Result Text
    let resultText = 'The system has a <span class="result-type">Unique Solution</span>:\n';
    for (let i = 0; i < M; i++) {
        resultText += `X${i+1} = ${formatNumber(X[i])}\n`;
    }
    const resultDiv = document.createElement('div');
    resultDiv.className = 'final-result';
    resultDiv.innerHTML = resultText;
    stepsContainer.appendChild(resultDiv);
}

// ====== Standard Helpers ======
function getEquationString(matrixRow, N) {
    let eq = '';
    let isFirstTerm = true;
    for (let j = 0; j < N; j++) {
        let coeff = matrixRow[j];
        if (Math.abs(coeff) < TOLERANCE) continue; 
        let coeffSign = coeff > 0 ? '+' : '-';
        let coeffAbs = Math.abs(coeff);
        let coeffText = formatNumber(coeffAbs);
        
        if (isFirstTerm) {
            coeffSign = coeff > 0 ? '' : '-';
            isFirstTerm = false;
        } else if (coeff > 0) {
            coeffSign = ' + ';
        } else {
            coeffSign = ' - ';
        }
        
        if (Math.abs(coeffAbs - 1) < TOLERANCE) {
            coeffText = '';
        }
        eq += `${coeffSign}${coeffText}x${j + 1}`;
    }
    eq += ` = ${formatNumber(matrixRow[N])}`;
    return eq.trim().replace(/^\+\s*/, '');
}

function backSubstitution(M, N, matrix, stepsContainer) {
    stepsContainer.innerHTML += '<h3>4. Back-Substitution:</h3>';
    const solution = new Array(N).fill(0);
    let pivotRows = [];
    for (let i = 0; i < M; i++) {
        let firstNonZeroCol = -1;
        for (let j = 0; j < N; j++) {
            if (Math.abs(matrix[i][j] - 1) < TOLERANCE) { 
                firstNonZeroCol = j;
                break;
            }
        }
        if (firstNonZeroCol !== -1) {
            pivotRows.push({row: i, col: firstNonZeroCol});
        }
    }
    for (let k = pivotRows.length - 1; k >= 0; k--) {
        const {row: i, col: j} = pivotRows[k]; 
        let rhs = matrix[i][N];
        let equation = getEquationString(matrix[i], N); 
        for (let l = j + 1; l < N; l++) {
            rhs -= matrix[i][l] * solution[l];
        }
        solution[j] = rhs;
        const stepDiv = document.createElement('div');
        stepDiv.className = 'back-sub-step';
        let equationDisplay = `<div class="equation-text">Equation R${i+1}: ${equation}</div>`;
        let finalSolution = `<div class="solution-line">X${j+1} = ${formatNumber(solution[j])}</div>`;
        stepDiv.innerHTML = equationDisplay + finalSolution;
        stepsContainer.appendChild(stepDiv);
    }
    return solution;
}

function interpretResults(M, N, matrix, stepsContainer, method) {
    let pivotCols = [];
    let resultText = '';
    let consistent = true;
    let finalSolution = new Array(N).fill(0); 
    
    for (let i = 0; i < M; i++) {
        let isZeroRow = true;
        let pivotJ = -1;
        for (let j = 0; j < N; j++) {
            if (Math.abs(matrix[i][j]) >= TOLERANCE) {
                isZeroRow = false;
                if (pivotJ === -1) pivotJ = j;
            }
        }
        if (pivotJ !== -1) pivotCols.push(pivotJ);
        if (isZeroRow && Math.abs(matrix[i][N]) >= TOLERANCE) {
            consistent = false;
            break;
        }
    }
    
    if (!consistent) {
        resultText = 'The system is <span class="result-type">Inconsistent</span>. No solution exists.';
    } else if (pivotCols.length < N) {
        resultText = `The system has <span class="result-type">Infinite Solutions</span>. (Number of free variables: ${N - pivotCols.length})`;
        let basicVars = pivotCols.map(j => `X${j + 1}`);
        let freeVars = Array.from({length: N}, (_, j) => j).filter(j => !pivotCols.includes(j)).map(j => `X${j + 1}`);
        let solutionString = `\nBasic Variables (Dependent): ${basicVars.join(', ')}\n`;
        solutionString += `Free Variables (Independent): ${freeVars.join(', ')} (Can be any value)\n`;
        resultText += solutionString;
    } else { 
        resultText = 'The system has a <span class="result-type">Unique Solution</span>:';
        if (method === 'gaussian') {
            finalSolution = backSubstitution(M, N, matrix, stepsContainer);
        } else { 
            for (let i = 0; i < N; i++) {
                let pivotRow = -1;
                for (let r = 0; r < M; r++) {
                    if (pivotCols.includes(i) && Math.abs(matrix[r][i] - 1) < TOLERANCE) {
                        pivotRow = r;
                        break;
                    }
                }
                if (pivotRow !== -1) finalSolution[i] = matrix[pivotRow][N];
            }
        }
        let solutionString = '\n';
        for (let i = 0; i < N; i++) {
            solutionString += `X${i+1} = ${formatNumber(finalSolution[i])}\n`;
        }
        resultText += solutionString;
    }
    const resultDiv = document.createElement('div');
    resultDiv.className = 'final-result';
    resultDiv.innerHTML = resultText;
    stepsContainer.appendChild(resultDiv);
}

function solveSystem() {
    const M = parseInt(document.getElementById('num-rows').value);
    const TotalCols = parseInt(document.getElementById('total-cols').value);
    const N = TotalCols - 1;
    const method = document.getElementById('elimination-method').value;

    const matrix = getMatrix(M, N); 
    const stepsContainer = document.getElementById('steps-container');
    stepsContainer.innerHTML = '';
    
    if (document.getElementById('A-0-0') === null) {
        stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Please click "Build Matrix" first.</p>';
        return;
    }

    if (!matrix) {
        stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Input Error. Please fill all fields with correct numbers.</p>';
        return;
    }

    // === Handle Inverse Method ===
    if (method === 'inverse') {
        if (M !== N) {
            stepsContainer.innerHTML += `<p class="final-result" style="color: red;">Error: Inverse Method requires a Square Matrix (Rows must equal Variables).<br>Current: ${M} Equations, ${N} Variables.</p>`;
            return;
        }
        stepsContainer.innerHTML += `<h3>Solving using Inverse Matrix Method (X = A⁻¹ * b):</h3>`;
        solveInverseMethod(M, N, matrix, stepsContainer);
        return;
    }
    
    stepsContainer.innerHTML += `<h3>Solving using ${method === 'gaussian' ? 'Gaussian Elimination' : 'Gauss-Jordan Elimination'}:</h3>`;
    stepsContainer.innerHTML += '<h3>1. Initial Augmented Matrix:</h3>';
    displayMatrix(matrix, M, N, stepsContainer);

    forwardElimination(M, TotalCols, matrix, stepsContainer);
    
    if (method === 'gauss-jordan') {
        stepsContainer.innerHTML += '<h3>2. Backward Elimination (To RREF):</h3>';
        backwardElimination(M, TotalCols, matrix, stepsContainer);
        stepsContainer.innerHTML += '<h3>3. Final Reduced Row Echelon Form (RREF):</h3>';
    } else {
        stepsContainer.innerHTML += '<h3>2. Final Row Echelon Form (REF):</h3>';
    }
    
    displayMatrix(matrix, M, N, stepsContainer);
    interpretResults(M, N, matrix, stepsContainer, method);
}

function clearAll() {
    document.getElementById('matrix-container').innerHTML = '';
    document.getElementById('steps-container').innerHTML = '';
    document.getElementById('num-rows').value = '3';
    document.getElementById('total-cols').value = '4'; 
    buildMatrix(); 
}

window.onload = buildMatrix;
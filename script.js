// Define tolerance for floating point comparisons
const TOLERANCE = 1e-9; 

// Helper function to format numbers: avoids .0000 for integers and formats decimals nicely.
function formatNumber(num) {
    if (Math.abs(num - Math.round(num)) < TOLERANCE) {
        return Math.round(num).toString();
    }
    return num.toFixed(4).replace(/\.?0+$/, "");
}

// Helper function to construct clean operation text
function getOperationText(row1, row2, factor, opType) {
    if (Math.abs(factor) < TOLERANCE && opType !== 'swap') { 
        return ''; 
    }
    
    if (opType === 'swap') {
        return `Swap R${row1} <-> R${row2}`;
    }

    if (opType === 'multiply') {
        return `R${row1} = ${formatNumber(factor)} * R${row1}`;
    }

    let factorAbs = Math.abs(factor);
    let sign = factor > 0 ? '-' : '+';
    
    let factorString = formatNumber(factorAbs);
    if (Math.abs(factorAbs - 1) < TOLERANCE) {
        factorString = ''; 
    }

    return `R${row1} = R${row1} ${sign} ${factorString} R${row2}`;
}

// Display the matrix in a step
function displayMatrix(matrix, M, N, container) {
    const grid = document.createElement('div');
    grid.className = 'matrix-display';
    grid.style.gridTemplateColumns = `repeat(${N + 1}, auto)`; 

    for (let i = 0; i < M; i++) {
        for (let j = 0; j < N + 1; j++) {
            const cell = document.createElement('div');
            if (matrix[i] && matrix[i][j] !== undefined) {
                cell.textContent = formatNumber(matrix[i][j]); 
            } else {
                 cell.textContent = '';
            }
            cell.className = 'matrix-cell';

            if (j === N) {
                cell.style.borderLeft = '3px solid #000';
            }
            grid.appendChild(cell);
        }
    }
    container.appendChild(grid);
}

// Display a full step (Operation text + Matrix)
function displayStep(operation, matrix, M, N) {
    if (operation === "") return;

    const stepsContainer = document.getElementById('steps-container');
    const stepDiv = document.createElement('div');
    stepDiv.className = 'step-item';

    const operationText = document.createElement('span');
    operationText.className = 'operation-text';
    operationText.textContent = operation;
    stepDiv.appendChild(operationText);

    displayMatrix(matrix, M, N, stepDiv);
    stepsContainer.appendChild(stepDiv);
}

// Build the input fields for the matrix
function buildMatrix() {
    const M = parseInt(document.getElementById('num-rows').value); // Rows (Equations)
    const TotalCols = parseInt(document.getElementById('total-cols').value); // Total Columns (A|b)
    const N = TotalCols - 1; // Coefficients Columns (Variables)
    
    const container = document.getElementById('matrix-container');
    container.innerHTML = ''; 

    if (M < 1 || N < 1 || M > 10 || TotalCols > 11) {
        alert("Please ensure rows (1-10) and total columns (2-11) are within limits, and total columns > 1.");
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

// Read the matrix values from the input fields
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

// Forward Elimination (to Row Echelon Form - REF)
function forwardElimination(M, N, matrix, stepsContainer) {
    let currentPivotRow = 0;

    for (let j = 0; j < N && currentPivotRow < M; j++) { 

        // 1. Partial Pivoting
        let pivotRow = currentPivotRow;
        for (let i = currentPivotRow + 1; i < M; i++) {
            if (Math.abs(matrix[i][j]) > Math.abs(matrix[pivotRow][j])) {
                pivotRow = i;
            }
        }
        
        if (Math.abs(matrix[pivotRow][j]) < TOLERANCE) { 
            continue; 
        }

        // Swap rows if necessary
        if (pivotRow !== currentPivotRow) {
            [matrix[currentPivotRow], matrix[pivotRow]] = [matrix[pivotRow], matrix[currentPivotRow]];
            displayStep(getOperationText(currentPivotRow + 1, pivotRow + 1, 0, 'swap'), matrix, M, N);
        }

        // 2. Zero out elements below the pivot
        for (let i = currentPivotRow + 1; i < M; i++) {
            let factor = matrix[i][j] / matrix[currentPivotRow][j];
            
            if (Math.abs(factor) < TOLERANCE) continue; 

            let operation = getOperationText(i + 1, currentPivotRow + 1, factor, 'subtract');

            for (let k = j; k < N + 1; k++) {
                matrix[i][k] = matrix[i][k] - factor * matrix[currentPivotRow][k];
            }
            displayStep(operation, matrix, M, N);
        }
        
        // 3. Set the current pivot element to 1 (Normalization - part of REF)
        let pivotValue = matrix[currentPivotRow][j];
        let factor = 1 / pivotValue;
        
        if (Math.abs(pivotValue - 1) > TOLERANCE) { 
            let operation = getOperationText(currentPivotRow + 1, 0, factor, 'multiply');
            
            for (let k = j; k < N + 1; k++) {
                matrix[currentPivotRow][k] = matrix[currentPivotRow][k] * factor;
            }
            displayStep(operation, matrix, M, N);
        }
        
        currentPivotRow++;
    }
    return currentPivotRow; // Return rank (number of pivot rows)
}

// Backward Elimination (for Gauss-Jordan)
function backwardElimination(M, N, matrix, stepsContainer) {
     let pivotCols = [];
     // Find pivot columns
    for (let i = 0; i < M; i++) {
        for (let j = 0; j < N; j++) {
            if (Math.abs(matrix[i][j] - 1) < TOLERANCE) {
                pivotCols.push({row: i, col: j});
                break;
            }
        }
    }

    // Zero out elements above the pivots
    for (const {row: i, col: j} of pivotCols) {
        for (let k = i - 1; k >= 0; k--) {
            let factor = matrix[k][j]; 
            
            if (Math.abs(factor) < TOLERANCE) continue;

            let operation = getOperationText(k + 1, i + 1, factor, 'subtract');

            for (let l = j; l < N + 1; l++) {
                matrix[k][l] = matrix[k][l] - factor * matrix[i][l];
            }
            displayStep(operation, matrix, M, N);
        }
    }
}

// Convert a matrix row to a human-readable equation string
function getEquationString(matrixRow, N) {
    let eq = '';
    let isFirstTerm = true;

    for (let j = 0; j < N; j++) {
        let coeff = matrixRow[j];
        if (Math.abs(coeff) < TOLERANCE) continue; // Skip zero coefficients

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

// Perform Back-Substitution and display steps
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

// Interpret and display final results
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
        if (isZeroRow && Math.abs(matrix[i][N]) >= TOLERANCE) consistent = false;
    }
    
    if (!consistent) {
        resultText = 'The system is <span class="result-type">Inconsistent</span>. No solution exists.';
    } else if (pivotCols.length < N) {
        resultText = `The system has <span class="result-type">Infinite Solutions</span>. (Number of free variables: ${N - pivotCols.length})`;
        let basicVars = pivotCols.map(j => `X${j + 1}`);
        let freeVars = Array.from({length: N}, (_, j) => j).filter(j => !pivotCols.includes(j)).map(j => `X${j + 1}`);
        resultText += `\nBasic Variables: ${basicVars.join(', ')}\nFree Variables: ${freeVars.join(', ')} (Can be any value)\n`;
    } else {
        resultText = 'The system has a <span class="result-type">Unique Solution</span>:';
        if (method === 'gaussian') {
            finalSolution = backSubstitution(M, N, matrix, stepsContainer);
        } else {
            for (let i = 0; i < N; i++) {
                let pivotRow = -1;
                for (let r = 0; r < M; r++) {
                    if (pivotCols.includes(i) && Math.abs(matrix[r][i] - 1) < TOLERANCE) pivotRow = r;
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

// Main Solver Function
function solveSystem() {
    const M = parseInt(document.getElementById('num-rows').value);
    const TotalCols = parseInt(document.getElementById('total-cols').value);
    const N = TotalCols - 1;
    const method = document.getElementById('elimination-method').value;

    const matrix = getMatrix(M, N); 
    const stepsContainer = document.getElementById('steps-container');
    stepsContainer.innerHTML = '';
    
    if (document.getElementById(`A-0-0`) === null) {
        stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Please click "Build Matrix" first.</p>';
        return;
    }

    if (!matrix) {
        stepsContainer.innerHTML += '<p class="final-result" style="color: red;">Input Error. Please fill all fields with correct numbers.</p>';
        return;
    }
    
    stepsContainer.innerHTML += `<h3>Solving using ${method === 'gaussian' ? 'Gaussian Elimination' : 'Gauss-Jordan Elimination'}:</h3>`;
    stepsContainer.innerHTML += '<h3>1. Initial Augmented Matrix:</h3>';
    displayMatrix(matrix, M, N, stepsContainer);

    forwardElimination(M, N, matrix, stepsContainer);
    
    if (method === 'gauss-jordan') {
        stepsContainer.innerHTML += '<h3>2. Backward Elimination (To RREF):</h3>';
        backwardElimination(M, N, matrix, stepsContainer);
        stepsContainer.innerHTML += '<h3>3. Final Reduced Row Echelon Form (RREF):</h3>';
    } else {
        stepsContainer.innerHTML += '<h3>2. Final Row Echelon Form (REF):</h3>';
    }
    
    displayMatrix(matrix, M, N, stepsContainer);

    interpretResults(M, N, matrix, stepsContainer, method);
}

// Clear all inputs and steps
function clearAll() {
    document.getElementById('matrix-container').innerHTML = '';
    document.getElementById('steps-container').innerHTML = '';
    document.getElementById('num-rows').value = '3';
    document.getElementById('total-cols').value = '4'; 
    buildMatrix(); 
}

// Run buildMatrix on page load
window.onload = buildMatrix;

[CmdletBinding()]
param(
    [switch]$Overwrite
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$clinicalPath = Join-Path $root 'clinical_ppv_quadas3_first_pass.csv'
$pairedPath = Join-Path $root 'paired_bc_bm_quadas3_first_pass.csv'
$summaryPath = Join-Path $root 'quadas3_first_pass_summary.csv'
$questionsPath = Join-Path $root 'quadas3_signaling_questions.csv'
$guidancePath = Join-Path $root 'quadas3_tailoring_guidance.md'
$xlsxPath = Join-Path $root 'quadas3_first_pass.xlsx'

foreach ($path in @($clinicalPath, $pairedPath, $summaryPath, $questionsPath, $guidancePath)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing input: $path" }
}
if ((Test-Path -LiteralPath $xlsxPath) -and -not $Overwrite) {
    throw "Workbook exists. Re-run with -Overwrite only after preserving any human-entered fields: $xlsxPath"
}

$clinical = @(Import-Csv -LiteralPath $clinicalPath)
$paired = @(Import-Csv -LiteralPath $pairedPath)
$summary = @(Import-Csv -LiteralPath $summaryPath)
$questions = @(Import-Csv -LiteralPath $questionsPath)
if ($clinical.Count -ne 26) { throw "Expected 26 clinical estimates; found $($clinical.Count)" }
if ($paired.Count -ne 9) { throw "Expected 9 paired estimates; found $($paired.Count)" }
if ($questions.Count -ne 20) { throw "Expected 20 signaling questions; found $($questions.Count)" }

$dataColumns = @($clinical[0].PSObject.Properties.Name)
if ((@($paired[0].PSObject.Properties.Name) -join '|') -ne ($dataColumns -join '|')) {
    throw 'Clinical and paired CSV schemas differ.'
}

$requiredColumns = @(
    'estimate_id','evidence_layer','study','location','year_or_period','estimate_positive_tested',
    'analysis_role','source_basis','participants_risk_of_bias','participants_applicability',
    'index_test_risk_of_bias','index_test_applicability','target_condition_risk_of_bias',
    'target_condition_applicability','analysis_risk_of_bias','overall_risk_of_bias',
    'overall_applicability','ai_first_pass','human_confirmation_required','jhk_decision',
    'jhk_review_date','jhk_notes','independent_reviewer','independent_review_date',
    'independent_decision','independent_notes','consensus_status'
)
foreach ($column in $requiredColumns) {
    if ($column -notin $dataColumns) { throw "Missing required column: $column" }
}

$questionDefinitions = @{}
foreach ($question in $questions) { $questionDefinitions[$question.question_id] = $question.signaling_question }
$definitions = foreach ($column in $dataColumns) {
    $definition = switch ($column) {
        'estimate_id' { 'Stable estimate-level identifier.' }
        'evidence_layer' { 'Clinical-definition PPV or paired blood/bone-marrow evidence layer.' }
        'study' { 'Standardized first-author surname and publication year; site suffix retained when one report contributes multiple estimates.' }
        'location' { 'Clinical cohort or paired-study location.' }
        'year_or_period' { 'Outbreak or study period.' }
        'estimate_positive_tested' { 'Human-readable positive/tested pair or paired 2x2 cells.' }
        'analysis_role' { 'Current approved modeling role before risk-of-bias sensitivity analyses.' }
        'source_basis' { 'Count provenance and source-audit note used for the first pass.' }
        'participants_risk_of_bias' { 'QUADAS-3 Participants-domain risk-of-bias judgment.' }
        'participants_rationale' { 'Reason for Participants-domain risk-of-bias judgment.' }
        'participants_applicability' { 'Participants-domain applicability concern.' }
        'participants_applicability_rationale' { 'Reason for Participants applicability judgment.' }
        'index_test_risk_of_bias' { 'QUADAS-3 Index-test-domain risk-of-bias judgment.' }
        'index_test_rationale' { 'Reason for Index-test-domain risk-of-bias judgment.' }
        'index_test_applicability' { 'Index-test-domain applicability concern.' }
        'index_test_applicability_rationale' { 'Reason for Index-test applicability judgment.' }
        'target_condition_risk_of_bias' { 'QUADAS-3 Target-condition-domain risk-of-bias judgment.' }
        'target_condition_rationale' { 'Reason for Target-condition-domain risk-of-bias judgment.' }
        'target_condition_applicability' { 'Target-condition-domain applicability concern.' }
        'target_condition_applicability_rationale' { 'Reason for Target-condition applicability judgment.' }
        'analysis_risk_of_bias' { 'QUADAS-3 Analysis-domain risk-of-bias judgment.' }
        'analysis_rationale' { 'Reason for Analysis-domain risk-of-bias judgment.' }
        'overall_risk_of_bias' { 'Overall estimate-level QUADAS-3 risk-of-bias judgment.' }
        'overall_applicability' { 'Overall estimate-level applicability concern.' }
        'overall_rationale' { 'Reason for the overall judgments.' }
        'ai_first_pass' { 'TRUE because Codex produced the initial structured judgment.' }
        'human_confirmation_required' { 'TRUE until JHK and an independent human reviewer complete appraisal.' }
        'jhk_decision' { 'JHK response to the AI first pass.' }
        'jhk_review_date' { 'JHK review date in YYYY-MM-DD format.' }
        'jhk_notes' { 'JHK corrections or supporting rationale.' }
        'independent_reviewer' { 'Name or initials of the independent human reviewer.' }
        'independent_review_date' { 'Independent review date in YYYY-MM-DD format.' }
        'independent_decision' { 'Independent reviewer response to the first-pass row.' }
        'independent_notes' { 'Independent reviewer corrections or rationale.' }
        'consensus_status' { 'Whether human judgments remain pending or have reached consensus.' }
        default {
            if ($questionDefinitions.ContainsKey($column)) { $questionDefinitions[$column] }
            else { "Field used in the QUADAS-3 first-pass record: $column" }
        }
    }
    $allowed = switch -Wildcard ($column) {
        'q?_?' { 'Y; PY; PN; N; NI' }
        '*risk_of_bias' { 'low; high; insufficient_information' }
        '*applicability' { 'low; high; insufficient_information' }
        'jhk_decision' { 'agree_with_first_pass; revise; needs_full_text; not_applicable' }
        'independent_decision' { 'agree_with_first_pass; revise; needs_full_text; not_applicable' }
        'consensus_status' { 'pending; agreed; resolved_after_discussion' }
        default { '' }
    }
    [pscustomobject]@{
        column_name = $column
        definition = $definition
        type_or_format = if ($column -match 'date$') { 'YYYY-MM-DD' } else { 'text' }
        allowed_values = $allowed
    }
}

$startRows = @(
    [pscustomobject]@{item='Status'; value='AI-prefilled QUADAS-3 first pass; JHK confirmation and independent human review required'},
    [pscustomobject]@{item='Tool version'; value='QUADAS-3 version 1.2 structure; appraisal created 2026-08-21'},
    [pscustomobject]@{item='Clinical estimates'; value='26 human-approved estimate-level rows'},
    [pscustomobject]@{item='Paired estimates'; value='9 JHK-approved paired blood/bone-marrow 2x2 rows'},
    [pscustomobject]@{item='Primary caution'; value='High risk does not automatically exclude an estimate; use judgments to guide sensitivity analysis and interpretation'},
    [pscustomobject]@{item='Clinical adaptation'; value='Clinical definition is the index test; latent true S. Typhi infection is the target condition; blood culture is specific but imperfectly sensitive'},
    [pscustomobject]@{item='Paired adaptation'; value='Blood culture is the index test; bone-marrow culture is informative but is not treated as an infallible gold standard'},
    [pscustomobject]@{item='JHK workflow'; value='Review each row, edit judgments if needed, complete jhk_decision/date/notes, then obtain an independent human review'},
    [pscustomobject]@{item='Editable first-pass workbook'; value='Do not rerun this builder after entering human decisions unless those fields have first been synchronized to the source CSVs'},
    [pscustomobject]@{item='Official source'; value='https://www.bristol.ac.uk/population-health-sciences/projects/quadas/quadas-3/quadas-3-tool/'},
    [pscustomobject]@{item='No model rerun'; value='No PPV, renewal, or ORI model was run while producing this appraisal'}
)

$guidanceRows = @(Get-Content -LiteralPath $guidancePath | ForEach-Object { [pscustomobject]@{guidance_line=$_} })

$excel = $null
$workbook = $null
try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $excel.DisplayAlerts = $false
    $workbook = $excel.Workbooks.Add()

    while ($workbook.Worksheets.Count -gt 1) { $workbook.Worksheets.Item(2).Delete() }
    $startSheet = $workbook.Worksheets.Item(1)
    $startSheet.Name = 'START_HERE'
    foreach ($name in @('clinical_ppv','paired_bc_bm','summary','signaling_questions','tailoring_guidance','data_dictionary')) {
        $sheet = $workbook.Worksheets.Add([Type]::Missing, $workbook.Worksheets.Item($workbook.Worksheets.Count))
        $sheet.Name = $name
    }

    function Write-Table($sheet, $rows, $columns, [int]$freezeColumns = 0) {
        $matrix = New-Object 'object[,]' ($rows.Count + 1), $columns.Count
        for ($c = 0; $c -lt $columns.Count; $c++) { $matrix[0, $c] = $columns[$c] }
        for ($r = 0; $r -lt $rows.Count; $r++) {
            for ($c = 0; $c -lt $columns.Count; $c++) {
                $value = $rows[$r].($columns[$c])
                $matrix[($r + 1), $c] = if ($null -eq $value) { '' } else { [string]$value }
            }
        }
        $lastRow = $rows.Count + 1
        $lastColumn = $columns.Count
        $range = $sheet.Range($sheet.Cells.Item(1,1), $sheet.Cells.Item($lastRow,$lastColumn))
        $range.NumberFormat = '@'
        $range.Value2 = $matrix
        $header = $sheet.Range($sheet.Cells.Item(1,1), $sheet.Cells.Item(1,$lastColumn))
        $header.Font.Bold = $true
        $header.Font.Color = 16777215
        $header.Interior.Color = 12611584
        $header.WrapText = $true
        $range.VerticalAlignment = -4160
        $range.WrapText = $true
        $range.AutoFilter() | Out-Null
        $sheet.Activate()
        $excel.ActiveWindow.SplitRow = 1
        $excel.ActiveWindow.SplitColumn = $freezeColumns
        $excel.ActiveWindow.FreezePanes = $true
        $range.Columns.AutoFit() | Out-Null
        for ($c = 1; $c -le $lastColumn; $c++) {
            if ($sheet.Columns.Item($c).ColumnWidth -gt 45) { $sheet.Columns.Item($c).ColumnWidth = 45 }
        }
        $sheet.Rows.Item(1).RowHeight = 42
    }

    Write-Table $startSheet $startRows @('item','value') 1
    Write-Table $workbook.Worksheets.Item('clinical_ppv') $clinical $dataColumns 8
    Write-Table $workbook.Worksheets.Item('paired_bc_bm') $paired $dataColumns 8
    Write-Table $workbook.Worksheets.Item('summary') $summary @('evidence_layer','judgment_type','domain','judgment','n_estimates') 3
    Write-Table $workbook.Worksheets.Item('signaling_questions') $questions @('question_id','domain','signaling_question','allowed_answers','answer_key') 2
    Write-Table $workbook.Worksheets.Item('tailoring_guidance') $guidanceRows @('guidance_line') 0
    Write-Table $workbook.Worksheets.Item('data_dictionary') @($definitions) @('column_name','definition','type_or_format','allowed_values') 1

    foreach ($sheetName in @('clinical_ppv','paired_bc_bm')) {
        $sheet = $workbook.Worksheets.Item($sheetName)
        $lastRow = $sheet.UsedRange.Rows.Count
        $headerMap = @{}
        for ($c = 1; $c -le $sheet.UsedRange.Columns.Count; $c++) { $headerMap[[string]$sheet.Cells.Item(1,$c).Text] = $c }

        foreach ($column in @('jhk_decision','independent_decision')) {
            $range = $sheet.Range($sheet.Cells.Item(2,$headerMap[$column]), $sheet.Cells.Item($lastRow,$headerMap[$column]))
            $range.Validation.Delete()
            $range.Validation.Add(3,1,1,'agree_with_first_pass,revise,needs_full_text,not_applicable')
            $range.Validation.InCellDropdown = $true
            $range.Interior.Color = 13434879
        }
        $consensusRange = $sheet.Range($sheet.Cells.Item(2,$headerMap['consensus_status']), $sheet.Cells.Item($lastRow,$headerMap['consensus_status']))
        $consensusRange.Validation.Delete()
        $consensusRange.Validation.Add(3,1,1,'pending,agreed,resolved_after_discussion')
        $consensusRange.Validation.InCellDropdown = $true

        foreach ($field in @('participants_risk_of_bias','index_test_risk_of_bias','target_condition_risk_of_bias','analysis_risk_of_bias','overall_risk_of_bias','participants_applicability','index_test_applicability','target_condition_applicability','overall_applicability')) {
            $columnIndex = $headerMap[$field]
            for ($r = 2; $r -le $lastRow; $r++) {
                $cell = $sheet.Cells.Item($r,$columnIndex)
                switch ([string]$cell.Text) {
                    'low' { $cell.Interior.Color = 13434828 }
                    'high' { $cell.Interior.Color = 13551615 }
                    'insufficient_information' { $cell.Interior.Color = 10092543 }
                }
            }
        }
    }

    $startSheet.Activate()
    $workbook.SaveAs($xlsxPath, 51)
    $workbook.Close($false)
    $workbook = $null
}
finally {
    if ($null -ne $workbook) { $workbook.Close($false) }
    if ($null -ne $excel) {
        $excel.Quit()
        [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($excel)
    }
    [gc]::Collect()
    [gc]::WaitForPendingFinalizers()
}

Write-Output "Created $xlsxPath"

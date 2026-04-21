/* ── API Response Types ───────────────────────────────────────── */

export interface MethodOption {
  key: string;
  name: string;
}

export interface StrategyInfo {
  number: number;
  name: string;
  description: string;
  available_tools: string[];
}

/* ── Job Types ───────────────────────────────────────────────── */

export interface AssemblyConfig {
  mode: 'assemble' | 'custom';
  bcell_csv_path: string;
  ctl_csv_path: string;
  htl_csv_path: string;
  bcell_count: number;
  ctl_count: number;
  htl_count: number;
  assembly_order: string;
  add_adjuvant: boolean;
  add_his_tag: boolean;
  custom_sequence: string;
  custom_fasta_path: string;
  run_sasa: boolean;
  sasa_csv_path: string;
}

export interface JobSubmitRequest {
  input_value: string;
  strategy: number;
  pathogen_type: string;
  mhci_method: string;
  mhcii_method: string;
  selected_tools?: string[];
  pre_predicted_fastas?: Record<string, string>;
  assembly_config?: AssemblyConfig;
}

export interface JobSubmitResponse {
  job_id: string;
  status: string;
}

export interface JobStatusResponse {
  job_id: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  progress_pct: number;
  current_tool: string;
  failed_tools: string[];
  error: string;
}

export interface OutputFile {
  relative_path: string;
  download_url: string;
  size_bytes: number;
}

export interface JobResultsResponse {
  job_id: string;
  status: string;
  results_dir: string;
  outputs: OutputFile[];
  failed_tools: string[];
}

/* ── Wizard State ────────────────────────────────────────────── */

export type InputMode = 'upload' | 'paste' | 'uniprot';

export type PathogenType = 'bacteria' | 'virus' | 'protozoa' | 'fungi';

export interface WizardState {
  /* Step 1 */
  inputMode: InputMode;
  inputValue: string;
  fileName: string;
  pathogenType: PathogenType;
  /* Step 2 */
  strategy: number | null;
  mhciMethod: string;
  mhciiMethod: string;
  selectedTools: string[];
  prePredictedFastas: Record<string, string>;
  assemblyConfig: AssemblyConfig;
}

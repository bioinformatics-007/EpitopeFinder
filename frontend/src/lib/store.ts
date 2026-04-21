import { create } from 'zustand';
import type { WizardState, InputMode, PathogenType, AssemblyConfig } from './types';

const DEFAULT_ASSEMBLY: AssemblyConfig = {
  mode: 'assemble',
  bcell_csv_path: '',
  ctl_csv_path: '',
  htl_csv_path: '',
  bcell_count: 0,
  ctl_count: 0,
  htl_count: 0,
  assembly_order: '1',
  add_adjuvant: false,
  add_his_tag: false,
  custom_sequence: '',
  custom_fasta_path: '',
  run_sasa: false,
  sasa_csv_path: '',
};

interface WizardStore extends WizardState {
  setInputMode: (mode: InputMode) => void;
  setInputValue: (value: string) => void;
  setFileName: (name: string) => void;
  setPathogenType: (type: PathogenType) => void;
  setStrategy: (strategy: number) => void;
  setMhciMethod: (method: string) => void;
  setMhciiMethod: (method: string) => void;
  toggleTool: (tool: string) => void;
  setPrePredictedFasta: (key: string, path: string) => void;
  updateAssembly: (updates: Partial<AssemblyConfig>) => void;
  reset: () => void;
}

const initialState: WizardState = {
  inputMode: 'upload',
  inputValue: '',
  fileName: '',
  pathogenType: 'bacteria',
  strategy: null,
  mhciMethod: 'f',
  mhciiMethod: 'nmel',
  selectedTools: [],
  prePredictedFastas: {},
  assemblyConfig: { ...DEFAULT_ASSEMBLY },
};

export const useWizardStore = create<WizardStore>((set) => ({
  ...initialState,

  setInputMode: (mode) => set({ inputMode: mode }),
  setInputValue: (value) => set({ inputValue: value }),
  setFileName: (name) => set({ fileName: name }),
  setPathogenType: (type) => set({ pathogenType: type }),
  setStrategy: (strategy) => set({ strategy }),
  setMhciMethod: (method) => set({ mhciMethod: method }),
  setMhciiMethod: (method) => set({ mhciiMethod: method }),

  toggleTool: (tool) =>
    set((state) => ({
      selectedTools: state.selectedTools.includes(tool)
        ? state.selectedTools.filter((t) => t !== tool)
        : [...state.selectedTools, tool],
    })),

  setPrePredictedFasta: (key, path) =>
    set((state) => ({
      prePredictedFastas: { ...state.prePredictedFastas, [key]: path },
    })),

  updateAssembly: (updates) =>
    set((state) => ({
      assemblyConfig: { ...state.assemblyConfig, ...updates },
    })),

  reset: () => set({ ...initialState, assemblyConfig: { ...DEFAULT_ASSEMBLY } }),
}));

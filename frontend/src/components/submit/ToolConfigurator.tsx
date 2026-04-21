"use client";

import React, { useEffect, useState } from 'react';
import { useWizardStore } from '@/lib/store';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/Card';
import { Tooltip } from '@/components/ui/Tooltip';
import { Select } from '@/components/ui/Select';
import { FileUpload } from '@/components/ui/FileUpload';
import { Input } from '@/components/ui/Input';
import { api } from '@/lib/api';
import type { MethodOption, StrategyInfo } from '@/lib/types';


export function ToolConfigurator() {
  const [mhciOptions, setMhciOptions] = useState<{value: string, label: string}[]>([]);
  const [mhciiOptions, setMhciiOptions] = useState<{value: string, label: string}[]>([]);
  const [strategyInfo, setStrategyInfo] = useState<StrategyInfo | null>(null);

  const { 
    strategy, 
    mhciMethod, setMhciMethod, 
    mhciiMethod, setMhciiMethod,
    selectedTools, toggleTool,
    prePredictedFastas, setPrePredictedFasta,
    assemblyConfig, updateAssembly
  } = useWizardStore();

  useEffect(() => {
    Promise.all([
      api.getMhciMethods(),
      api.getMhciiMethods(),
      api.getStrategies()
    ]).then(([mhcI, mhcII, strats]) => {
      setMhciOptions(mhcI.map(m => ({ value: m.key, label: `${m.key} - ${m.name}` })));
      setMhciiOptions(mhcII.map(m => ({ value: m.key, label: `${m.key.toUpperCase()} - ${m.name}` })));
      if (strategy) {
        setStrategyInfo(strats.find(s => s.number === strategy) || null);
      }
    });
  }, [strategy]);

  if (!strategy) return null;

  const showMhcSelects = [1, 5, 2].includes(strategy); // Strategy 2 might need it if bcell/mhc selected
  const showToolCheckboxes = [1, 2, 3].includes(strategy);

  return (
    <div className="space-y-6 pb-24">
      {/* ── Methods Configuration ─────────────────────── */}
      {showMhcSelects && (
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center">
              Prediction Algorithms
              <Tooltip text="Algorithms used to predict binding affinity to MHC-I (Cytotoxic T Cells) and MHC-II (Helper T Cells) alleles." />
            </CardTitle>
            <CardDescription>Select the algorithms used for epitope prediction.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-6 md:grid-cols-2">
            <Select
              label="MHC-I Method"
              tooltip="Major Histocompatibility Complex Class I: essential for presenting endogenous peptides to CD8+ T cells."
              options={mhciOptions}
              value={mhciMethod}
              onChange={(e) => setMhciMethod(e.target.value)}
            />
            <Select
              label="MHC-II Method"
              tooltip="Major Histocompatibility Complex Class II: essential for presenting exogenous peptides to CD4+ T cells."
              options={mhciiOptions}
              value={mhciiMethod}
              onChange={(e) => setMhciiMethod(e.target.value)}
            />
          </CardContent>
        </Card>
      )}

      {/* ── Tool Selection (Strats 1, 2, 3) ────────────── */}
      {showToolCheckboxes && strategyInfo && (
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center">
              Target Prediction Tools
              <Tooltip text="B-cell tools predict linear epitopes. T-cell tools like NetCTL predict proteasomal cleavage and TAP transport." />
            </CardTitle>
            <CardDescription>Select the specific biological analysis tools to execute in this job.</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
              {strategyInfo.available_tools.map((tool) => {
                const isSelected = selectedTools.includes(tool) || selectedTools.length === 0; 
                return (
                  <button
                    key={tool}
                    onClick={() => toggleTool(tool)}
                    className={`flex items-center justify-between rounded-sm border-2 p-4 text-left ${
                      isSelected
                        ? 'border-amber-600 bg-amber-50/50 text-amber-900 '
                        : 'border-stone-100 bg-stone-50 text-stone-500 hover:bg-white hover:text-stone-800'
                    }`}
                  >
                    <span className="text-sm font-bold ">{tool}</span>
                    <div className={`flex h-6 w-6 items-center justify-center border-2 ${
                      isSelected ? 'border-amber-600 bg-amber-600 ' : 'border-stone-200 bg-white'
                    }`}>
                      
                    </div>
                  </button>
                );
              })}
            </div>
          </CardContent>
        </Card>
      )}

      {/* ── Strategy 5 specific (Pre-predicted files) ── */}
      {strategy === 5 && (
        <Card>
          <CardHeader>
            <CardTitle>Biological Data Ingress</CardTitle>
            <CardDescription>Ingest previously generated epitope results for secondary filtering and assembly.</CardDescription>
          </CardHeader>
          <CardContent className="space-y-6">
            <FileUpload
              label="B-Cell Epitopes FASTA"
              selectedFile={prePredictedFastas.bcell ? new File([''], prePredictedFastas.bcell) : null}
              onFileSelect={(f) => setPrePredictedFasta('bcell', f?.name || '')}
            />
            <FileUpload
              label="MHC-I Epitopes FASTA"
              selectedFile={prePredictedFastas.mhci ? new File([''], prePredictedFastas.mhci) : null}
              onFileSelect={(f) => setPrePredictedFasta('mhci', f?.name || '')}
            />
            <FileUpload
              label="MHC-II Epitopes FASTA"
              selectedFile={prePredictedFastas.mhcii ? new File([''], prePredictedFastas.mhcii) : null}
              onFileSelect={(f) => setPrePredictedFasta('mhcii', f?.name || '')}
            />
          </CardContent>
        </Card>
      )}

      {/* ── Strategy 6 specific (Assembly Config) ────── */}
      {strategy === 6 && (
        <Card>
          <CardHeader>
            <CardTitle>Vaccine Architecture Design</CardTitle>
            <CardDescription>Define the parameters for multi-epitope construct assembly and structural validation.</CardDescription>
          </CardHeader>
          <CardContent className="space-y-8">
            <div className="flex gap-4 p-1.5 rounded-sm bg-stone-100 border border-stone-200 w-fit">
               {(['assemble', 'custom'] as const).map(mode => (
                 <button 
                   key={mode}
                   onClick={() => updateAssembly({ mode })}
                   className={`px-6 py-2 rounded-sm text-sm font-bold ${
                     assemblyConfig.mode === mode 
                      ? 'bg-white text-amber-700 ' 
                      : 'text-stone-500 hover:text-stone-700'
                   }`}
                 >
                   <span className="capitalize">{mode} sequence</span>
                 </button>
               ))}
            </div>

            {assemblyConfig.mode === 'assemble' && (
              <div className="space-y-8">
                <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                  <div className="space-y-2">

                    <Input 
                      type="number"
                      value={assemblyConfig.bcell_count} 
                      onChange={e => updateAssembly({ bcell_count: parseInt(e.target.value) || 0 })} 
                    />
                  </div>
                  <div className="space-y-2">

                    <Input 
                      type="number"
                      value={assemblyConfig.ctl_count} 
                      onChange={e => updateAssembly({ ctl_count: parseInt(e.target.value) || 0 })} 
                    />
                  </div>
                  <div className="space-y-2">

                    <Input 
                      type="number" 
                      value={assemblyConfig.htl_count} 
                      onChange={e => updateAssembly({ htl_count: parseInt(e.target.value) || 0 })} 
                    />
                  </div>
                </div>
                
                <div className="space-y-2 pt-6 border-t border-stone-100">
                  <label className="text-sm font-bold text-stone-900 flex items-center">
                    Assembly Order
                    <Tooltip text="Provide comma separated order, e.g. adju_bcell_ctl_htl. Default is 1." />
                  </label>
                  <Input 
                    placeholder="e.g. 1" 
                    value={assemblyConfig.assembly_order}
                    onChange={e => updateAssembly({ assembly_order: e.target.value })}
                  />
                </div>
              </div>
            )}

            {assemblyConfig.mode === 'custom' && (
              <div className="space-y-2">
                 <label className="text-sm font-bold text-stone-900">Custom Construct Sequence</label>
                 <textarea
                  className="w-full rounded-sm border border-stone-200 bg-stone-50 px-4 py-3 text-sm text-stone-900 font-mono min-h-[160px] focus:bg-white focus:focus:outline-none "
                  placeholder="Enter protein sequence (amino acids only)..."
                  value={assemblyConfig.custom_sequence}
                  onChange={(e) => updateAssembly({ custom_sequence: e.target.value })}
                />
              </div>
            )}

            <div className="grid grid-cols-1 md:grid-cols-3 gap-6 pt-8 border-t border-stone-100">
              <label className={`flex items-center gap-3 p-4 rounded-sm border-2 cursor-pointer ${
                assemblyConfig.add_adjuvant ? 'border-amber-600 bg-amber-50/50' : 'border-stone-100 bg-stone-50'
              }`}>
                <input 
                  type="checkbox" 
                  className="h-5 w-5 rounded border-stone-300 text-amber-600 focus:"
                  checked={assemblyConfig.add_adjuvant}
                  onChange={e => updateAssembly({ add_adjuvant: e.target.checked })}
                />
                <span className="text-sm font-bold text-stone-900">Add Adjuvant</span>
              </label>
              
              <label className={`flex items-center gap-3 p-4 rounded-sm border-2 cursor-pointer ${
                assemblyConfig.add_his_tag ? 'border-amber-600 bg-amber-50/50' : 'border-stone-100 bg-stone-50'
              }`}>
                <input 
                  type="checkbox" 
                  className="h-5 w-5 rounded border-stone-300 text-amber-600 focus:"
                  checked={assemblyConfig.add_his_tag}
                  onChange={e => updateAssembly({ add_his_tag: e.target.checked })}
                />
                <span className="text-sm font-bold text-stone-900">Add 6xHis Tag</span>
              </label>

              <label className={`flex items-center gap-3 p-4 rounded-sm border-2 cursor-pointer ${
                assemblyConfig.run_sasa ? 'border-amber-600 bg-amber-50/50' : 'border-stone-100 bg-stone-50'
              }`}>
                <input 
                  type="checkbox" 
                  className="h-5 w-5 rounded border-stone-300 text-amber-600 focus:"
                  checked={assemblyConfig.run_sasa}
                  onChange={e => updateAssembly({ run_sasa: e.target.checked })}
                />
                <span className="text-sm font-bold text-stone-900">SASA Analysis</span>
              </label>
            </div>
          </CardContent>
        </Card>
      )}

    </div>
  );
}

"use client";

import React, { useState } from 'react';
import { useRouter } from 'next/navigation';
import { useWizardStore } from '@/lib/store';
import { api } from '@/lib/api';
import { Card, CardHeader, CardTitle, CardContent } from '@/components/ui/Card';
import { Button } from '@/components/ui/Button';
import { Badge } from '@/components/ui/Badge';


export function ReviewSubmit() {
  const router = useRouter();
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const state = useWizardStore();

  const handleSubmit = async () => {
    setIsSubmitting(true);
    setError(null);
    try {
      let jobId = '';
      
      if (state.inputMode === 'upload' && state.inputValue) {
        // We stored the raw text in inputValue, but API expects a file. 
        // We reconstruct a File object to send via FormData
        const formData = new FormData();
        const file = new File([state.inputValue], state.fileName || 'sequence.fasta', { type: 'text/plain' });
        formData.append('file', file);
        formData.append('strategy', String(state.strategy));
        formData.append('pathogen_type', state.pathogenType);
        formData.append('mhci_method', state.mhciMethod);
        formData.append('mhcii_method', state.mhciiMethod);
        if (state.selectedTools.length > 0) formData.append('selected_tools', JSON.stringify(state.selectedTools));
        if (state.strategy === 5) formData.append('pre_predicted_fastas', JSON.stringify(state.prePredictedFastas));
        if (state.strategy === 6) formData.append('assembly_config', JSON.stringify(state.assemblyConfig));

        const res = await api.submitWithFile(formData);
        jobId = res.job_id;
      } else {
        const payload = {
          input_value: state.inputValue,
          strategy: state.strategy!,
          pathogen_type: state.pathogenType,
          mhci_method: state.mhciMethod,
          mhcii_method: state.mhciiMethod,
          selected_tools: state.selectedTools.length > 0 ? state.selectedTools : undefined,
          pre_predicted_fastas: state.strategy === 5 ? state.prePredictedFastas : undefined,
          assembly_config: state.strategy === 6 ? state.assemblyConfig : undefined,
        };
        const res = await api.submitJob(payload);
        jobId = res.job_id;
      }

      router.push(`/results/${jobId}`);
    } catch (err: any) {
      console.error(err);
      setError(err.message || 'Failed to submit job.');
      setIsSubmitting(false);
    }
  };

  return (
    <div className="space-y-8 fade-in slide-in-">
      <Card className="border-amber-200 bg-amber-50/20 overflow-visible relative">
        <div className="absolute -top-4 -right-4 bg-amber-600 text-white p-3 rounded-sm ">
          
        </div>

        <CardHeader className="pb-6 border-b border-amber-100/50">
          <div className="flex flex-col space-y-1">
            <CardTitle className="text-amber-900 text-2xl">Configuration Finalized</CardTitle>
            <p className="text-sm font-medium text-amber-700/70">Review your biological assessment parameters before execution.</p>
          </div>
        </CardHeader>
        <CardContent className="pt-8 space-y-8">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="space-y-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="h-1.5 w-1.5 rounded-sm bg-amber-600" />

              </div>
              <div className="rounded-sm bg-white border border-stone-100 p-6 space-y-4 group ">
                <dl className="space-y-4">
                  <div className="flex justify-between items-center">
                    <dt className="text-sm font-bold text-stone-500">Input Mode</dt>

                  </div>
                  <div className="flex justify-between items-center pt-3 border-t border-stone-50">
                    <dt className="text-sm font-bold text-stone-500">Target Type</dt>

                  </div>
                  {state.inputMode === 'upload' && (
                    <div className="flex justify-between items-center pt-3 border-t border-stone-50">
                      <dt className="text-sm font-bold text-stone-500">Source File</dt>
                      <dd className="text-sm font-bold text-stone-900 truncate max-w-[180px]">{state.fileName}</dd>
                    </div>
                  )}
                  {state.inputMode === 'uniprot' && (
                    <div className="flex justify-between items-center pt-3 border-t border-stone-50">
                      <dt className="text-sm font-bold text-stone-500">UniProt ID</dt>
                      <dd className="text-sm font-bold text-amber-600 font-mono text-lg">{state.inputValue}</dd>
                    </div>
                  )}
                </dl>
              </div>
            </div>

            <div className="space-y-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="h-1.5 w-1.5 rounded-sm bg-amber-600" />

              </div>
              <div className="rounded-sm bg-white border border-stone-100 p-6 space-y-4 group ">
                <dl className="space-y-4">
                  <div className="flex justify-between items-center">
                    <dt className="text-sm font-bold text-stone-500">Strategy</dt>
                    <dd className="flex items-center gap-2">
                      <span className="h-7 w-7 flex items-center justify-center bg-amber-900 text-white text-xs font-bold">{state.strategy}</span>
                    </dd>
                  </div>
                  <div className="flex justify-between items-center pt-3 border-t border-stone-50">
                    <dt className="text-sm font-bold text-stone-500">MHC-I Protocol</dt>

                  </div>
                  <div className="flex justify-between items-center pt-3 border-t border-stone-50">
                    <dt className="text-sm font-bold text-stone-500">MHC-II Protocol</dt>

                  </div>
                </dl>
              </div>
            </div>

            {(state.selectedTools.length > 0 || state.strategy === 6) && (
              <div className="lg:col-span-2 space-y-4">
                <div className="flex items-center gap-2 mb-2">
                  <div className="h-1.5 w-1.5 rounded-sm bg-amber-600" />

                </div>
                <div className="rounded-sm bg-white border border-stone-100 p-6 group flex flex-wrap gap-4">
                  {state.selectedTools.map(t => (
                    <Badge key={t} variant="success" className="px-4 py-1.5 rounded-sm border-amber-100 bg-amber-50/50 text-amber-800 text-[11px] font-bold ">{t}</Badge>
                  ))}
                  {state.strategy === 6 && (
                    <>
                      <Badge variant="default" className="px-4 py-1.5 rounded-sm border-stone-200 bg-stone-50 text-stone-600 text-[11px] font-bold ">Mode: {state.assemblyConfig.mode}</Badge>
                      {state.assemblyConfig.add_adjuvant && <Badge variant="success" className="px-4 py-1.5 rounded-sm border-amber-100 bg-amber-500 text-white text-[11px] font-bold ">+ Adjuvant</Badge>}
                      {state.assemblyConfig.add_his_tag && <Badge variant="success" className="px-4 py-1.5 rounded-sm border-amber-100 bg-amber-500 text-white text-[11px] font-bold ">+ 6xHis Tag</Badge>}
                      {state.assemblyConfig.run_sasa && <Badge variant="success" className="px-4 py-1.5 rounded-sm border-amber-100 bg-amber-500 text-white text-[11px] font-bold ">+ SASA</Badge>}
                    </>
                  )}
                </div>
              </div>
            )}
          </div>

          {error && (
            <div className="flex items-center gap-4 rounded-sm border-2 border-red-100 bg-red-50 p-6 text-red-700 ">
              
              <div className="flex flex-col">

                <p className="text-sm font-medium">{error}</p>
              </div>
            </div>
          )}

          <div className="flex justify-end pt-6">
            <Button 
              size="lg" 
              onClick={handleSubmit} 
              isLoading={isSubmitting} 
              className="w-full h-16 rounded-sm text-lg font-bold ] .02] active:scale-[0.98] bg-amber-600 hover:bg-amber-700 text-white border border-[#999]"
            >
              Initialize Assessment Pipeline
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}

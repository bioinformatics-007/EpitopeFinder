"use client";

import React, { useRef } from 'react';
import { useWizardStore } from '@/lib/store';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/Card';
import { Tabs } from '@/components/ui/Tabs';
import { Input } from '@/components/ui/Input';
import { FileUpload } from '@/components/ui/FileUpload';


export function SequenceInput() {
  const { inputMode, inputValue, pathogenType, setInputMode, setInputValue, setFileName, setPathogenType } = useWizardStore();
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileUpload = (file: File | null) => {
    if (!file) {
      setInputValue('');
      setFileName('');
      return;
    }
    setFileName(file.name);
    const reader = new FileReader();
    reader.onload = (e) => {
      const text = e.target?.result as string;
      setInputValue(text);
    };
    reader.readAsText(file);
  };

  const tabs = [
    { id: 'upload', label: 'Upload FASTA', },
    { id: 'paste', label: 'Paste Sequence', },
    { id: 'uniprot', label: 'UniProt ID', },
  ];

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>Sequence Input</CardTitle>
          <CardDescription>Provide the pathogen sequence to begin analysis.</CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          <Tabs
            tabs={tabs}
            activeTab={inputMode}
            onChange={(id) => {
              setInputMode(id as any);
              setInputValue('');
              setFileName('');
            }}
          />

          <div className="min-h-[200px]">
            {inputMode === 'upload' && (
              <div className="space-y-4">
                <FileUpload
                  label="Pathogen File"
                  accept=".fasta,.fa,.txt"
                  selectedFile={inputValue ? new File([inputValue], useWizardStore.getState().fileName || 'sequence.fasta') : null}
                  onFileSelect={handleFileUpload}
                />
                {!inputValue && (
                  <button
                    onClick={() => {
                      setFileName('test_monkeypox.fasta');
                      setInputValue('>sp|A0A7H0DNC4|PG153_MONPV Envelop protein OPG153 OS=Monkeypox virus OX=10244 GN=OPG153 PE=3 SV=1\nMANIINLWNGIVPMVQDVNVASITAFKSMIDETWDKKIEANTCISRKHRNIIHEVIRDFM\nKAYPKMDENRKSPLGAPMQWLTQYYILKNEYHKTMLAYDDGSLNTKFKTLNIYMITNVGQ\nYILYIVFCIISGKNHDGTPYIYDSEITSNDKNLINDRIKYACKQILHGQLTMALRIRNKF\nMFIGSPMYLWFNVNGSHVYHEIYDRNVGFHNKEIGRLLYAFMYYLSISGRFLNDLALLKF\nTYLGESWTFSLSVPEYILYGLGYSVFDTIEKFSNDAILVYIRTNNRNGYDYAEFNKKGIV\nKVTEDKPDNDKRIHAIRLINYSSDVQHIHFGFRNTLIIDNECTNIQSSAENATDIGHYQD\nSKINIEDDDIIDDDDDDDDDDDDDDYNPKPTPIPDPHPRPPFPRHDYHKRPKLLPVEEPD\nPVKKDADRIRLDNHILNTLDHNLNSIGHYCCDTVAVDRLEHHIETLGQYTVILARKINMQ\nTLLFPWPLPTVHQHAIDGSIPPHGRSTIL');
                      setPathogenType('virus');
                    }}
                    className="flex items-center gap-2 text-xs font-bold text-amber-600 hover:text-amber-700 "
                  >
                     Use Test Data (Monkeypox Virus)
                  </button>
                )}
              </div>
            )}

            {inputMode === 'paste' && (
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <label className="text-sm font-bold text-stone-900">Raw Sequence data</label>
                  <button
                    onClick={() => {
                      setInputValue('>sp|A0A7H0DNC4|PG153_MONPV Envelop protein OPG153 OS=Monkeypox virus OX=10244 GN=OPG153 PE=3 SV=1\nMANIINLWNGIVPMVQDVNVASITAFKSMIDETWDKKIEANTCISRKHRNIIHEVIRDFM\nKAYPKMDENRKSPLGAPMQWLTQYYILKNEYHKTMLAYDDGSLNTKFKTLNIYMITNVGQ\nYILYIVFCIISGKNHDGTPYIYDSEITSNDKNLINDRIKYACKQILHGQLTMALRIRNKF\nMFIGSPMYLWFNVNGSHVYHEIYDRNVGFHNKEIGRLLYAFMYYLSISGRFLNDLALLKF\nTYLGESWTFSLSVPEYILYGLGYSVFDTIEKFSNDAILVYIRTNNRNGYDYAEFNKKGIV\nKVTEDKPDNDKRIHAIRLINYSSDVQHIHFGFRNTLIIDNECTNIQSSAENATDIGHYQD\nSKINIEDDDIIDDDDDDDDDDDDDDYNPKPTPIPDPHPRPPFPRHDYHKRPKLLPVEEPD\nPVKKDADRIRLDNHILNTLDHNLNSIGHYCCDTVAVDRLEHHIETLGQYTVILARKINMQ\nTLLFPWPLPTVHQHAIDGSIPPHGRSTIL');
                      setPathogenType('virus');
                    }}
                    className="text-xs font-bold text-amber-600 hover:text-amber-700 flex items-center gap-1"
                  >
                     Use Test Sequence
                  </button>
                </div>
                <textarea
                  className="w-full flex-1 rounded-sm border border-stone-200 bg-stone-50 px-4 py-3 text-sm text-stone-900 placeholder:text-stone-400 focus-visible:outline-none focus-visible:focus-visible:font-mono min-h-[200px] focus:bg-white"
                  placeholder=">Sequence_1\nMKVLY..."
                  value={inputValue}
                  onChange={(e) => setInputValue(e.target.value)}
                  spellCheck={false}
                />
              </div>
            )}

            {inputMode === 'uniprot' && (
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <label className="text-sm font-bold text-stone-900">UniProt Accession ID</label>
                  <button
                    onClick={() => {
                      setInputValue('A0A7H0DNC4');
                      setPathogenType('virus');
                    }}
                    className="text-xs font-bold text-amber-600 hover:text-amber-700 flex items-center gap-1"
                  >
                     Use Test ID (Monkeypox)
                  </button>
                </div>
                <Input
                  placeholder="e.g. P0DTC2"
                  value={inputValue}
                  onChange={(e) => setInputValue(e.target.value.toUpperCase())}
                  className="font-mono text-xl h-14 border-stone-200 focus:rounded-sm"
                />
                <p className="text-xs text-stone-500 font-medium">
                  The system will automatically fetch the FASTA sequence from UniProt.
                </p>
              </div>
            )}
          </div>

          <div className="space-y-4 pt-8 border-t border-stone-200">

            <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
              {(['bacteria', 'virus', 'protozoa', 'fungi'] as const).map((type) => (
                <button
                  key={type}
                  onClick={() => setPathogenType(type)}
                  className={`flex items-center justify-center rounded-sm border-2 px-4 py-4 text-sm font-bold ${
                    pathogenType === type
                      ? 'border-amber-600 bg-amber-50 text-amber-700 '
                      : 'border-stone-100 bg-stone-50 text-stone-500 hover:bg-white hover:text-stone-800'
                  }`}
                >
                  <span className="capitalize">{type}</span>
                </button>
              ))}
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}

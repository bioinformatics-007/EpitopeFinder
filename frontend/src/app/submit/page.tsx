'use client';

import React, { useState } from 'react';
import { useWizardStore } from '@/lib/store';
import { Button } from '@/components/ui/Button';
import { SequenceInput } from '@/components/submit/SequenceInput';
import { StrategySelector } from '@/components/submit/StrategySelector';
import { ToolConfigurator } from '@/components/submit/ToolConfigurator';
import { ReviewSubmit } from '@/components/submit/ReviewSubmit';


const STEPS = [
  { id: 1, title: 'Input Data' },
  { id: 2, title: 'Strategy' },
  { id: 3, title: 'Configuration' },
  { id: 4, title: 'Review' },
];

export default function SubmitWizardPage() {
  const [currentStep, setCurrentStep] = useState(1);
  const state = useWizardStore();

  // Validation logic to prevent moving forward if invalid
  const canGoNext = () => {
    if (currentStep === 1) {
      if (!state.inputValue) return false;
      return true;
    }
    if (currentStep === 2) {
      if (!state.strategy) return false;
      return true;
    }
    if (currentStep === 3) {
      // Configuration step is mostly pre-filled with defaults. Just return true.
      return true;
    }
    return false;
  };

  const handleNext = () => {
    if (canGoNext() && currentStep < 4) {
      setCurrentStep(s => s + 1);
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
  };

  const handlePrev = () => {
    if (currentStep > 1) {
      setCurrentStep(s => s - 1);
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
  };

  return (
    <div className="text-sm text-black leading-normal">
      {/* ── Step Indicator ─────────────────────────────────── */}
      <div className="mb-6 border border-[#999]">
        <table className="w-full border-collapse text-sm">
          <tbody>
            <tr>
              {STEPS.map((step) => {
                const isCompleted = currentStep > step.id;
                const isCurrent = currentStep === step.id;
                return (
                  <td
                    key={step.id}
                    className={`border border-[#999] px-3 py-2 text-center font-bold ${
                      isCurrent
                        ? 'bg-[#3A332D] text-white'
                        : isCompleted
                          ? 'bg-[#F4F4EE] text-[#3A332D]'
                          : 'bg-white text-stone-400'
                    }`}
                  >
                    {isCompleted ? '✓ ' : ''}{step.title}
                  </td>
                );
              })}
            </tr>
          </tbody>
        </table>
      </div>

      {/* ── Step Content ───────────────────────────────────── */}
      <div className="min-h-[300px] border border-[#999] p-4 bg-white">
        {currentStep === 1 && <SequenceInput />}
        {currentStep === 2 && <StrategySelector />}
        {currentStep === 3 && <ToolConfigurator />}
        {currentStep === 4 && <ReviewSubmit />}
      </div>

      {/* ── Navigation Buttons ─────────────────────────────── */}
      <div className="mt-4 flex justify-between items-center">
        {currentStep > 1 ? (
          <button onClick={handlePrev} className="px-4 py-2 border border-[#999] bg-[#F4F4EE] text-sm font-bold text-[#3A332D] cursor-pointer hover:bg-[#e5e5d0]">
            ← Previous Step
          </button>
        ) : <div />}
        
        {currentStep < 4 ? (
          <button  
            onClick={handleNext} 
            disabled={!canGoNext()} 
            className="px-6 py-2 bg-[#3A332D] text-white text-sm font-bold border border-[#999] cursor-pointer hover:bg-[#635B53] disabled:opacity-50 disabled:cursor-not-allowed"
          >
            Next Step →
          </button>
        ) : <div />}
      </div>
    </div>
  );
}

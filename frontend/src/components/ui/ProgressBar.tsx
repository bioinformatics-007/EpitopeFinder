import React from 'react';

interface ProgressBarProps {
  value: number; // 0 to 100
  className?: string;
  showLabel?: boolean;
}

export function ProgressBar({ value, className = '', showLabel = false }: ProgressBarProps) {
  // Clamp value between 0 and 100
  const safeValue = Math.min(100, Math.max(0, value));

  return (
    <div className={`w-full ${className}`}>
      {showLabel && (
        <div className="mb-2 flex justify-between text-xs font-bold ">
          <span className="text-stone-900">Execution Progress</span>
          <span className="text-amber-600 bg-amber-50 px-2 py-0.5 border border-amber-100">{Math.round(safeValue)}%</span>
        </div>
      )}
      <div className="h-3 w-full overflow-hidden rounded-sm bg-stone-100 border border-stone-200 p-0.5 ">
        <div
          className="h-full ease-out 0_0_10px_rgba(16,185,129,0.3)] relative overflow-hidden"
          style={{ width: `${safeValue}%` }}
        >
          <div className="absolute inset-0 bg-[linear-gradient(45deg,rgba(255,255,255,0.15)_25%,transparent_25%,transparent_50%,rgba(255,255,255,0.15)_50%,rgba(255,255,255,0.15)_75%,transparent_75%,transparent)] bg-[length:1rem_1rem] stripe_1s_linear_infinite]" />
        </div>
      </div>
    </div>
  );
}

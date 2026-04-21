import React from 'react';
import { Tooltip } from '@/components/ui/Tooltip';

export interface SelectProps extends React.SelectHTMLAttributes<HTMLSelectElement> {
  options: { value: string; label: string }[];
  label?: string;
  error?: string;
  tooltip?: string;
}

export const Select = React.forwardRef<HTMLSelectElement, SelectProps>(
  ({ className = '', options, label, error, tooltip, ...props }, ref) => {
    return (
      <div className="flex flex-col space-y-1.5 w-full">
        {label && (
          <label className="text-sm font-bold text-stone-900 flex items-center mb-1">
            {label}
            {tooltip && <Tooltip text={tooltip} />}
          </label>
        )}
        <div className="relative">
          <select
            className={`flex h-12 w-full appearance-none rounded-sm border border-stone-200 bg-stone-50 px-4 py-2 text-sm text-stone-900 focus:bg-white focus:focus:focus:border-amber-600 outline-none disabled:cursor-not-allowed disabled:opacity-50 pr-10 ${
              error ? 'border-red-500 focus:' : ''
            } ${className}`}
            ref={ref}
            {...props}
          >
            {options.map((opt) => (
              <option key={opt.value} value={opt.value} className="bg-white text-stone-900">
                {opt.label}
              </option>
            ))}
          </select>
          <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center pr-3">
            <svg className="h-5 w-5 text-stone-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
            </svg>
          </div>
        </div>
        {error && <p className="text-sm text-red-400">{error}</p>}
      </div>
    );
  }
);
Select.displayName = 'Select';

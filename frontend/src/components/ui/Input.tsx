import React from 'react';

export interface InputProps extends React.InputHTMLAttributes<HTMLInputElement> {
  error?: string;
}

export const Input = React.forwardRef<HTMLInputElement, InputProps>(
  ({ className = '', type, error, ...props }, ref) => {
    return (
      <div className="flex flex-col space-y-1.5 w-full">
        <input
          type={type}
          className={`flex h-12 w-full rounded-sm border border-stone-200 bg-stone-50 px-4 py-2 text-sm text-stone-900 placeholder:text-stone-400 focus:bg-white focus:focus:focus:border-amber-600 outline-none disabled:cursor-not-allowed disabled:opacity-50 ${
            error ? 'border-red-500 focus:' : ''
          } ${className}`}
          ref={ref}
          {...props}
        />
        {error && <p className="text-sm text-red-400">{error}</p>}
      </div>
    );
  }
);
Input.displayName = 'Input';

import React from 'react';

type BadgeVariant = 'default' | 'success' | 'warning' | 'error' | 'outline';

interface BadgeProps extends React.HTMLAttributes<HTMLDivElement> {
  variant?: BadgeVariant;
}

export const Badge = React.forwardRef<HTMLDivElement, BadgeProps>(
  ({ className = '', variant = 'default', children, ...props }, ref) => {
    const baseStyles =
      'inline-flex items-center rounded-sm px-2.5 py-0.5 text-[10px] font-bold focus:outline-none focus:focus:focus:';

    const variants: Record<BadgeVariant, string> = {
      default: 'bg-stone-100 text-stone-600 border border-stone-200',
      success: 'bg-amber-50 text-amber-700 border border-amber-200 ',
      warning: 'bg-amber-50 text-amber-700 border border-amber-200',
      error: 'bg-red-50 text-red-700 border border-red-200',
      outline: 'text-stone-500 border border-stone-300',
    };

    return (
      <div ref={ref} className={`${baseStyles} ${variants[variant]} ${className}`} {...props}>
        {children}
      </div>
    );
  }
);

Badge.displayName = 'Badge';

import React from 'react';


export function Tooltip({ text }: { text: string }) {
  return (
    <div className="group relative inline-flex items-center ml-2">
      
      <div className="pointer-events-none absolute bottom-full left-1/2 mb-2 hidden -transtone-x-1/2 w-max max-w-xs rounded-sm bg-stone-900 px-3 py-2 text-[11px] font-medium text-white ">
        {text}
        <div className="absolute left-1/2 top-full -mt-[1px] -transtone-x-1/2 border-[6px] border-transparent border-t-stone-900" />
      </div>
    </div>
  );
}

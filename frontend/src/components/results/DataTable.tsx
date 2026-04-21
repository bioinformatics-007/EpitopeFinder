"use client";

import React, { useState, useEffect } from 'react';


interface DataTableProps {
  csvUrl: string;
  filename: string;
}

export function DataTable({ csvUrl, filename }: DataTableProps) {
  const [data, setData] = useState<Record<string, string>[]>([]);
  const [columns, setColumns] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);
  const [search, setSearch] = useState('');
  const [sortCol, setSortCol] = useState<string | null>(null);
  const [sortDesc, setSortDesc] = useState(false);

  useEffect(() => {
    fetch(csvUrl)
      .then(res => res.text())
      .then(csv => {
        // Simple CSV parser
        const lines = csv.split('\n').map(l => l.trim()).filter(l => l.length > 0);
        if (lines.length > 0) {
          const headers = lines[0].split(',').map(h => h.replace(/^"|"$/g, '').trim());
          const rows = lines.slice(1).map(line => {
            // Very naive split, doesn't handle commas inside quotes perfectly, 
            // but usually sufficient for bioinformatics tabular outputs unless there are complex descriptions.
            const values = line.split(',').map(v => v.replace(/^"|"$/g, '').trim());
            const rowObj: Record<string, string> = {};
            headers.forEach((h, i) => { rowObj[h] = values[i] || ''; });
            return rowObj;
          });
          setColumns(headers);
          setData(rows);
        }
      })
      .catch(console.error)
      .finally(() => setLoading(false));
  }, [csvUrl]);

  if (loading) {
    return <div className="h-64 w-full skeleton " />;
  }

  if (data.length === 0) {
    return (
    <div className="flex h-32 items-center justify-center rounded-sm border-2 border-dashed border-stone-200 bg-stone-50/50 m-4">

    </div>
  );
}

// Filtering
const filteredData = data.filter(row => 
  Object.values(row).some(val => val.toLowerCase().includes(search.toLowerCase()))
);

// Sorting
const sortedData = [...filteredData].sort((a, b) => {
  if (!sortCol) return 0;
  const valA = parseFloat(a[sortCol]) || a[sortCol];
  const valB = parseFloat(b[sortCol]) || b[sortCol];
  
  if (valA < valB) return sortDesc ? 1 : -1;
  if (valA > valB) return sortDesc ? -1 : 1;
  return 0;
});

const handleSort = (col: string) => {
  if (sortCol === col) {
    setSortDesc(!sortDesc);
  } else {
    setSortCol(col);
    setSortDesc(false);
  }
};

return (
  <div className="flex flex-col flex-1 h-full max-h-full">
    <div className="flex items-center justify-between px-6 py-4 border-b border-stone-100 bg-white">
      <div className="relative w-80 group">
        
        <input
          type="text"
          placeholder="Filter clinical records..."
          value={search}
          onChange={e => setSearch(e.target.value)}
          className="h-11 w-full rounded-sm border border-stone-200 bg-stone-50 pl-11 pr-4 text-sm text-stone-900 placeholder:text-stone-400 focus:bg-white focus:border-amber-600 focus:outline-none focus:focus:font-medium"
        />
      </div>
      <a 
        href={csvUrl} 
        download={filename}
        className="flex items-center gap-2 rounded-sm bg-stone-100 px-4 py-2 text-xs font-bold text-stone-700 hover:bg-amber-600 hover:text-white hover:hover:"
      >
         Export CSV
      </a>
    </div>

    <div className="overflow-auto bg-white flex-1 relative custom-scrollbar">
      <table className="w-full text-left text-sm text-stone-700 whitespace-nowrap min-w-full">
        <thead className="sticky top-0 bg-stone-50 text-[10px] font-bold ] text-stone-400 ">
          <tr className="border-b border-stone-100">
            {columns.map(col => (
              <th 
                key={col} 
                className="px-6 py-4 font-bold cursor-pointer hover:text-amber-700 group"
                onClick={() => handleSort(col)}
              >
                <div className="flex items-center gap-2">
                  {col}
                  <div className="flex flex-col justify-center gap-0.5 opacity-0 ">
                     
                     
                  </div>
                </div>
              </th>
            ))}
          </tr>
        </thead>
        <tbody className="divide-y divide-stone-100 font-mono text-xs">
          {sortedData.map((row, i) => (
            <tr key={i} className="data-row group hover:bg-amber-50/30">
              {columns.map(col => (
                <td key={col} className="px-6 py-3.5 text-stone-600 border-r border-stone-50/50 last:border border-[#999]">
                  {row[col]}
                </td>
              ))}
            </tr>
          ))}
          {sortedData.length === 0 && (
            <tr>
              <td colSpan={columns.length} className="px-6 py-24 text-center">
                <div className="flex flex-col items-center gap-4 opacity-50">
                   

                </div>
              </td>
            </tr>
          )}
        </tbody>
      </table>
    </div>
  </div>
);
}

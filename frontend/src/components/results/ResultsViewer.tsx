"use client";

import React, { useState } from 'react';
import type { OutputFile } from '@/lib/types';
import { Card, CardHeader, CardTitle, CardContent } from '@/components/ui/Card';
import { Tabs } from '@/components/ui/Tabs';
import { DataTable } from '@/components/results/DataTable';

import { Button } from '@/components/ui/Button';

interface ResultsViewerProps {
  outputs: OutputFile[];
  jobId: string;
}

export function ResultsViewer({ outputs, jobId }: ResultsViewerProps) {
  const [activeTab, setActiveTab] = useState(outputs[0]?.relative_path || '');

  // Filter for actual data files (CSVs) we want to tab-view
  // You might want to view txt files as raw text blocks, but we'll stick to CSV in DataTable for now
  const viewableFiles = outputs.filter(o => o.relative_path.endsWith('.csv') || o.relative_path.endsWith('.txt'));

  return (
    <Card className="flex flex-col h-[800px] border-stone-200 overflow-hidden mt-12 fade-in slide-in-bg-white ">
      <CardHeader className="border-b border-stone-100 bg-stone-50/50 px-8 py-6 shrink-0 relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-full opacity-50" />
        <div className="flex items-center justify-between relative ">
          <div className="flex items-center gap-4">
            <div className="rounded-sm bg-amber-600 p-3 text-white ">
              
            </div>
            <div>
              <CardTitle className="text-2xl font-bold text-stone-900 ">Analytical Insights</CardTitle>

            </div>
          </div>
          
          <Button variant="outline" size="sm" className="rounded-sm font-bold border-stone-200 hover:bg-amber-50 hover:text-amber-700 ">
             Download All (ZIP)
          </Button>
        </div>
      </CardHeader>
      
      <div className="flex flex-1 overflow-hidden">
        {/* Sidebar Nav */}
        <div className="w-80 shrink-0 overflow-y-auto border-r border-stone-100 bg-stone-50 px-4 py-8 space-y-6">
          <div>

            <div className="space-y-1.5">
              {viewableFiles.map(file => {
                const isActive = activeTab === file.relative_path;
                const fileName = file.relative_path.split('/').pop() || '';
                const pathParts = file.relative_path.split('/');
                const category = pathParts.length > 1 ? pathParts[0] : 'General';
                
                return (
                  <button
                    key={file.relative_path}
                    onClick={() => setActiveTab(file.relative_path)}
                    className={`flex w-full items-center gap-4 rounded-sm px-4 py-4 text-sm group ${
                      isActive
                        ? 'bg-white text-amber-900 font-bold border border-stone-100 transtone-x-1'
                        : 'text-stone-500 hover:bg-white hover:text-stone-800 hover:'
                    }`}
                  >
                    <div className={`p-2 ${
                      isActive ? 'bg-amber-600 text-white' : 'bg-stone-200 text-stone-400 '
                    }`}>
                      
                    </div>
                    <div className="flex flex-col items-start truncate text-left">
                      <span className="truncate w-full font-bold ">{fileName}</span>
                      <span className={`text-[10px] font-bold ${isActive ? 'text-amber-600/70' : 'text-stone-400 '}`}>
                        {category}
                      </span>
                    </div>
                  </button>
                );
              })}
            </div>
          </div>
        </div>

        {/* Main Content Area */}
        <div className="flex-1 bg-white flex flex-col min-w-0">
          <div className="flex-1 p-6 overflow-hidden">
            {activeTab ? (
              <div className="h-full border border-stone-100 rounded-sm overflow-hidden bg-stone-50/30">
                <DataTable 
                  csvUrl={outputs.find(o => o.relative_path === activeTab)!.download_url} 
                  filename={activeTab.split('/').pop() || 'data.csv'} 
                />
              </div>
            ) : (
              <div className="flex flex-col h-full items-center justify-center text-stone-300 space-y-4">
                

              </div>
            )}
          </div>
        </div>
      </div>
    </Card>
  );
}

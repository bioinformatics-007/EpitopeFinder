"use client";

import React, { useEffect, useState } from 'react';
import { api } from '@/lib/api';
import type { JobStatusResponse } from '@/lib/types';
import { Card, CardContent } from '@/components/ui/Card';
import { ProgressBar } from '@/components/ui/ProgressBar';
import { Badge } from '@/components/ui/Badge';


interface JobProgressProps {
  jobId: string;
  onComplete: () => void;
  onError: (msg: string) => void;
}

export function JobProgress({ jobId, onComplete, onError }: JobProgressProps) {
  const [status, setStatus] = useState<JobStatusResponse | null>(null);

  useEffect(() => {
    let interval: NodeJS.Timeout;

    const poll = async () => {
      try {
        const data = await api.getJobStatus(jobId);
        setStatus(data);

        if (data.status === 'completed') {
          onComplete();
          return true; // stop polling
        } else if (data.status === 'failed') {
          onError(data.error || 'Pipeline failed during execution.');
          return true; // stop polling
        }
      } catch (err) {
        console.error('Polling error', err);
      }
      return false; // keep polling
    };

    // Initial poll
    poll().then((stop) => {
      if (!stop) {
        interval = setInterval(async () => {
          const shouldStop = await poll();
          if (shouldStop) clearInterval(interval);
        }, 3000);
      }
    });

    return () => clearInterval(interval);
  }, [jobId, onComplete, onError]);

  if (!status) {
    return (
      <Card className="border-stone-200 bg-stone-50/50 border-dashed">
        <CardContent className="flex flex-col items-center justify-center p-16">
          <div className="relative mb-6">
            <div className="absolute inset-0 bg-amber-400/20 rounded-sm" />
          </div>
          <h3 className="text-xl font-bold text-stone-900 mb-2">Synchronizing Protocol</h3>
          <p className="text-sm font-medium text-stone-500">Establishing secure handshake with the execution environment...</p>
        </CardContent>
      </Card>
    );
  }

  const isRunning = status.status === 'running' || status.status === 'pending';
  const isCompleted = status.status === 'completed';
  const isFailed = status.status === 'failed';

  return (
    <Card className={`border-stone-200 ${isRunning ? 'border-amber-300' : isFailed ? 'border-red-200' : ''}`}>
      {/* Active Assessment banner — sits ABOVE the card, not overlapping the badge */}
      {isRunning && (
        <div className="w-full flex justify-center -mt-px">
          <div className="px-5 py-1 bg-amber-600 text-white text-[10px] font-bold tracking-widest uppercase rounded-t-none rounded-b-sm border border-amber-700">
            Active Assessment
          </div>
        </div>
      )}

      <CardContent className="p-10">
        <div className="flex flex-col md:flex-row items-center md:items-start justify-between gap-8 mb-10">
          <div className="flex flex-col md:flex-row items-center gap-6 text-center md:text-left">
            {/* Status label box */}
            <div className={`p-4 rounded-sm whitespace-nowrap font-bold text-sm ${
              isCompleted ? 'bg-amber-600 text-white' :
              isFailed    ? 'bg-red-600 text-white' :
              'bg-white border-2 border-amber-600 text-amber-600'
            }`}>
              {isCompleted ? '[ Complete ]' : isFailed ? '[ Failed ]' : '[ Processing ]'}
            </div>

            <div>
              {/* Title + status badge — on the same row, no absolute positioning */}
              <div className="flex flex-wrap items-center gap-3 mb-2">
                <h2 className="text-2xl font-bold text-stone-900">
                  Pipeline Execution
                </h2>
                <Badge variant={
                  isCompleted ? 'success' : isFailed ? 'error' : 'default'
                } className="h-7 px-4">
                  {isRunning && <span className="h-2 w-2 rounded-sm bg-amber-400 mr-2 inline-block" />}
                  <span className="font-bold">{status.status}</span>
                </Badge>
              </div>
              <p className="text-xs text-stone-400 font-bold flex items-center gap-2 justify-center md:justify-start">
                Reference ID: <span className="text-stone-600 font-mono text-[13px] normal-case">{jobId}</span>
              </p>
            </div>
          </div>
        </div>

        <div className="space-y-6">
          <ProgressBar
            value={status.progress_pct}
            showLabel
            className={isFailed ? 'grayscale opacity-50' : ''}
          />

          {/* Current tool box — always shown (not just when running) */}
          <div className="rounded-sm bg-stone-50 p-6 border border-stone-100">
            <div className="flex items-center gap-4">
              {isRunning && (
                <div className="relative h-10 w-10 flex items-center justify-center bg-white rounded-sm border border-stone-100 shrink-0">
                  <span className="absolute inline-flex h-4 w-4 rounded-sm bg-amber-400 opacity-75 animate-ping" />
                  <div className="relative h-2.5 w-2.5 rounded-sm bg-amber-600" />
                </div>
              )}
              <span className="font-mono text-sm font-bold text-amber-800 bg-amber-50/50 px-4 py-2 border border-amber-100/50">
                {status.current_tool || 'Initializing secure environment...'}
              </span>
            </div>
          </div>

          {/* Failed tools warning */}
          {status.failed_tools?.length > 0 && (
            <div className="rounded-sm bg-amber-50/30 border border-amber-200/50 p-6">
              <div className="flex flex-wrap gap-2">
                {status.failed_tools.map(tool => (
                  <Badge key={tool} variant="warning" className="bg-amber-100/50 text-amber-700 border-amber-200">{tool} alert</Badge>
                ))}
              </div>
            </div>
          )}

          {/* Error box — clearly separated, full width */}
          {isFailed && (
            <div className="rounded-sm bg-red-50 border-2 border-red-100 p-6 mt-2">
              <p className="text-xs font-bold uppercase tracking-widest text-red-400 mb-2">Execution Error</p>
              <p className="text-sm font-medium text-red-700 font-mono break-all">
                {status.error || 'The biological assessment was intercepted by an unhandled exception.'}
              </p>
            </div>
          )}
        </div>
      </CardContent>
    </Card>
  );
}

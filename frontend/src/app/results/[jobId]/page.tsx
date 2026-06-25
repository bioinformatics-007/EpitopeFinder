'use client';

import React, { useEffect, useState } from 'react';
import { api } from '@/lib/api';
import type { JobResultsResponse } from '@/lib/types';
import { JobProgress } from '@/components/results/JobProgress';
import { ResultsViewer } from '@/components/results/ResultsViewer';


export default function ResultsPage({ params }: { params: { jobId: string } }) {
  const { jobId } = params;
  const [results, setResults] = useState<JobResultsResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  const fetchResults = async () => {
    try {
      const data = await api.getJobResults(jobId);
      setResults(data);
    } catch (err: any) {
      setError(err.message || 'Failed to fetch job results.');
    }
  };

  return (
    <div className="mx-auto max-w-6xl px-4 py-12 sm:px-6 lg:px-8">

      {/* Top-level page error banner (network/fetch failures) */}
      {error && !results && (
        <div className="mb-8 rounded-sm border-2 border-red-200 bg-red-50 p-5">
          <p className="text-xs font-bold uppercase tracking-widest text-red-400 mb-1">Pipeline Error</p>
          <p className="text-sm text-red-700 font-mono break-all">{error}</p>
        </div>
      )}

      {/* Progress tracking — shown while running or failed; hidden once results arrive */}
      {!results && !error && (
        <JobProgress
          jobId={jobId}
          onComplete={fetchResults}
          onError={setError}
        />
      )}

      {/* Results dashboard — shown once job is completed and outputs fetched */}
      {results && results.outputs && (
        <div>
          <ResultsViewer outputs={results.outputs} jobId={jobId} />
        </div>
      )}
    </div>
  );
}

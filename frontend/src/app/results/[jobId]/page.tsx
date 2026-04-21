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

  // On mount, see if we already have results (maybe job is already done?)
  // We can just rely on JobProgress to poll first, but it's simpler to let JobProgress do its thing.
  // Actually, we'll wait for JobProgress to say 'complete' and then fetch the actual outputs.

  return (
    <div className="mx-auto max-w-6xl px-4 py-12 sm:px-6 lg:px-8">
      {error && (
        <div className="mb-8 flex items-center gap-2 border border-red-500/20 bg-red-500/10 p-4 text-red-400">
          
          <p>{error}</p>
        </div>
      )}

      {/* Progress tracking - hidden once results are fetched */}
      {!results && (
        <JobProgress 
          jobId={jobId} 
          onComplete={fetchResults}
          onError={setError}
        />
      )}

      {/* Results dashboard - shown once job is completed */}
      {results && results.outputs && (
        <div className="fade-in slide-in-">
          <ResultsViewer outputs={results.outputs} jobId={jobId} />
        </div>
      )}
    </div>
  );
}

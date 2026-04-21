import type {
  MethodOption,
  StrategyInfo,
  JobSubmitRequest,
  JobSubmitResponse,
  JobStatusResponse,
  JobResultsResponse,
} from './types';

const API_BASE = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

/* ── Helpers ─────────────────────────────────────────────────── */

async function get<T>(path: string): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`);
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`GET ${path} failed (${res.status}): ${body}`);
  }
  return res.json();
}

async function post<T>(path: string, data: unknown): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`POST ${path} failed (${res.status}): ${body}`);
  }
  return res.json();
}

async function postForm<T>(path: string, formData: FormData): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    method: 'POST',
    body: formData,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`POST ${path} failed (${res.status}): ${body}`);
  }
  return res.json();
}

/* ── API Client ──────────────────────────────────────────────── */

export const api = {
  /* Config */
  getMhciMethods: () => get<MethodOption[]>('/api/config/mhci-methods'),
  getMhciiMethods: () => get<MethodOption[]>('/api/config/mhcii-methods'),
  getStrategies: () => get<StrategyInfo[]>('/api/config/strategies'),

  /* Jobs */
  submitJob: (data: JobSubmitRequest) =>
    post<JobSubmitResponse>('/api/jobs/submit', data),

  submitWithFile: (formData: FormData) =>
    postForm<JobSubmitResponse>('/api/jobs/submit-with-file', formData),

  getJobStatus: (jobId: string) =>
    get<JobStatusResponse>(`/api/jobs/${jobId}/status`),

  getJobResults: (jobId: string) =>
    get<JobResultsResponse>(`/api/jobs/${jobId}/results`),

  getDownloadUrl: (jobId: string, filePath: string) =>
    `${API_BASE}/api/jobs/${jobId}/results/${filePath}`,

  /* Health */
  health: () => get<{ status: string; version: string }>('/api/health'),
};

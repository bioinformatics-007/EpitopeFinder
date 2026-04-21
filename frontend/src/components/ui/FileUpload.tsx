import React, { useCallback, useRef } from 'react';


interface FileUploadProps {
  accept?: string;
  onFileSelect: (file: File | null) => void;
  selectedFile: File | null;
  error?: string;
  label: string;
}

export function FileUpload({ accept = '.fasta,.fa,.txt', onFileSelect, selectedFile, error, label }: FileUploadProps) {
  const inputRef = useRef<HTMLInputElement>(null);

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      onFileSelect(e.target.files[0]);
    }
  };

  const clearFile = (e: React.MouseEvent) => {
    e.stopPropagation();
    onFileSelect(null);
    if (inputRef.current) inputRef.current.value = '';
  };

  return (
    <div className="flex flex-col space-y-1.5 w-full">
      <label className="text-sm font-medium text-stone-900">{label}</label>
      
      {!selectedFile ? (
        <div
          onClick={() => inputRef.current?.click()}
          className={`flex cursor-pointer flex-col items-center justify-center border-2 border-dashed bg-stone-100/50 px-6 py-8 hover:bg-stone-200/50 ${
            error ? 'border-red-500/50' : 'border-stone-300'
          }`}
        >
          
          <p className="text-sm font-medium text-stone-800">Click to upload or drag and drop</p>
          <p className="text-xs text-stone-500 mt-1">Accepts {accept.split(',').join(', ')}</p>
        </div>
      ) : (
        <div className="flex items-center justify-between border border-stone-300 bg-white px-4 py-3">
          <div className="flex items-center space-x-3 overflow-hidden">
            <div className="flex h-10 w-10 shrink-0 items-center justify-center rounded-md bg-amber-600/10">
              
            </div>
            <div className="flex flex-col truncate">
              <span className="truncate text-sm font-medium text-stone-900">{selectedFile.name}</span>
              <span className="text-xs text-stone-500">{(selectedFile.size / 1024).toFixed(1)} KB</span>
            </div>
          </div>
          <button
            onClick={clearFile}
            className="ml-4 rounded-sm p-1 text-stone-600 hover:bg-stone-200 hover:text-red-400"
            aria-label="Remove file"
          >
            
          </button>
        </div>
      )}
      
      <input
        type="file"
        ref={inputRef}
        onChange={handleFileChange}
        accept={accept}
        className="hidden"
      />
      
      {error && <p className="text-sm text-red-400">{error}</p>}
    </div>
  );
}

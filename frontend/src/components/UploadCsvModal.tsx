import React, { useRef, useState } from 'react';

interface Environment {
  environment_id: string;
  title: string;
}

interface UploadCsvModalProps {
  open: boolean;
  onClose: () => void;
  environments: Environment[];
  defaultEnvId?: string;
  onUploadSuccess?: () => void;
  refreshEnvFiles?: () => void;
}

const UploadCsvModal: React.FC<UploadCsvModalProps> = ({ open, onClose, environments, defaultEnvId, onUploadSuccess, refreshEnvFiles }) => {
  const [uploadEnvId, setUploadEnvId] = useState<string>(defaultEnvId || (environments[0]?.environment_id ?? ''));
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const [uploading, setUploading] = useState(false);
  const [dragActive, setDragActive] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  React.useEffect(() => {
    setUploadEnvId(defaultEnvId || (environments[0]?.environment_id ?? ''));
  }, [defaultEnvId, environments]);

  if (!open) return null;

  const handleUploadConfirm = async () => {
    if (!uploadFile) {
      setUploadError('Please select a file.');
      return;
    }
    setUploading(true);
    setUploadError(null);
    const formData = new FormData();
    formData.append('file', uploadFile);
    formData.append('environment_id', uploadEnvId);
    try {
      const res = await fetch(`/api/environments/${uploadEnvId}/files`, {
        method: 'POST',
        credentials: 'include',
        cache: 'no-store',
        body: formData,
      });
      if (!res.ok) throw new Error('Upload failed');
      setUploadFile(null);
      setUploadError(null);
      if (onUploadSuccess) onUploadSuccess();
      if (refreshEnvFiles) refreshEnvFiles();
      onClose();
    } catch (err) {
      setUploadError('Failed to upload file.');
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black bg-opacity-60">
      <div className="bg-gray-900 rounded-lg p-6 w-full max-w-md relative">
        <button
          className="absolute top-2 right-2 text-gray-400 hover:text-white"
          onClick={onClose}
          disabled={uploading}
        >
          ×
        </button>
        <h2 className="text-lg font-bold text-white mb-4">Upload CSV File</h2>
        <label className="block text-white mb-2">Select Environment</label>
        <select
          className="w-full p-2 rounded border mb-4 text-white bg-gray-800"
          value={uploadEnvId}
          onChange={e => setUploadEnvId(e.target.value)}
          disabled={uploading}
        >
          {(environments ?? []).map(env => (
            <option key={env.environment_id} value={env.environment_id}>{env.title}</option>
          ))}
        </select>
        {/* Drag and drop zone */}
        <div
          onDragOver={e => { e.preventDefault(); setDragActive(true); }}
          onDragLeave={e => { e.preventDefault(); setDragActive(false); }}
          onDrop={e => {
            e.preventDefault();
            setDragActive(false);
            if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
              const file = e.dataTransfer.files[0];
              if (file.type === 'text/csv' || file.name.endsWith('.csv')) {
                setUploadFile(file);
                setUploadError(null);
              } else {
                setUploadError('Please drop a valid CSV file.');
              }
            }
          }}
          onClick={() => fileInputRef.current?.click()}
          className={`w-full mb-4 p-4 border-2 rounded-lg transition-colors duration-200 ${dragActive ? 'border-green-400 bg-green-900/20' : 'border-dashed border-gray-600 bg-gray-800/40'}`}
          style={{ textAlign: 'center', color: dragActive ? '#22c55e' : '#fff', cursor: 'pointer' }}
        >
          {uploadFile ? (
            <span>Selected file: <span className="font-semibold text-green-400">{uploadFile.name}</span></span>
          ) : (
            <span>Drag and drop your CSV file here, or <span className="underline">click to browse</span>.</span>
          )}
        </div>
        <input
          type="file"
          accept=".csv"
          style={{ display: 'none' }}
          ref={fileInputRef}
          onChange={e => setUploadFile(e.target.files?.[0] || null)}
          disabled={uploading}
        />
        {uploadError && <div className="text-red-500 mb-2">{uploadError}</div>}
        <button
          className="w-full bg-green-700 text-white py-2 rounded hover:bg-green-800 disabled:opacity-50"
          onClick={handleUploadConfirm}
          disabled={uploading}
        >
          {uploading ? 'Uploading...' : 'Upload'}
        </button>
      </div>
    </div>
  );
};

export default UploadCsvModal; 
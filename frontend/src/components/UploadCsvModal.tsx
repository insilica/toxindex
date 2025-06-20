import React, { useState } from 'react';
import { FaTimes } from 'react-icons/fa';

interface Environment {
  environment_id: string;
  title: string;
}

interface UploadCsvModalProps {
  open: boolean;
  onClose: () => void;
  environments?: Environment[];
  defaultEnvId?: string;
  onUploadSuccess?: () => void;
  refreshEnvFiles?: (envId: string) => Promise<void>;
}

const UploadCsvModal: React.FC<UploadCsvModalProps> = ({
  open,
  onClose,
  environments = [],
  defaultEnvId,
  onUploadSuccess,
  refreshEnvFiles
}) => {
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [selectedEnv, setSelectedEnv] = useState(defaultEnvId || '');
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  if (!open) return null;

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      if (!file.name.toLowerCase().endsWith('.csv')) {
        setError('Please select a CSV file');
        return;
      }
      setSelectedFile(file);
      setError(null);
    }
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!selectedFile || !selectedEnv) {
      setError('Please select both a file and an environment');
      return;
    }

    setUploading(true);
    setError(null);

    const formData = new FormData();
    formData.append('file', selectedFile);

    try {
      const response = await fetch(`/api/environments/${selectedEnv}/files`, {
        method: 'POST',
        credentials: 'include',
        body: formData,
      });

      const data = await response.json();

      if (data.success) {
        if (onUploadSuccess) {
          onUploadSuccess();
        }
        if (refreshEnvFiles) {
          await refreshEnvFiles(selectedEnv);
        }
        onClose();
      } else {
        setError(data.error || 'Failed to upload file');
      }
    } catch (err) {
      setError('Failed to upload file');
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-gray-900 rounded-lg p-6 max-w-md w-full mx-4">
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-xl font-semibold text-white">Upload CSV File</h2>
          <button
            onClick={onClose}
            className="text-gray-400 hover:text-white transition-colors"
          >
            <FaTimes />
          </button>
        </div>

        <form onSubmit={handleSubmit}>
          <div className="mb-4">
            <label className="block text-gray-300 mb-2">Environment</label>
            <select
              value={selectedEnv}
              onChange={(e) => setSelectedEnv(e.target.value)}
              className="w-full bg-gray-800 text-white rounded px-3 py-2 focus:outline-none focus:ring-2 focus:ring-green-500"
              required
            >
              <option value="">Select an environment</option>
              {environments.map((env) => (
                <option key={env.environment_id} value={env.environment_id}>
                  {env.title}
                </option>
              ))}
            </select>
          </div>

          <div className="mb-4">
            <label className="block text-gray-300 mb-2">CSV File</label>
            <input
              type="file"
              accept=".csv"
              onChange={handleFileChange}
              className="w-full bg-gray-800 text-white rounded px-3 py-2 focus:outline-none focus:ring-2 focus:ring-green-500"
              required
            />
          </div>

          {error && (
            <div className="mb-4 text-red-500">{error}</div>
          )}

          <div className="flex justify-end">
            <button
              type="button"
              onClick={onClose}
              className="mr-2 px-4 py-2 text-gray-300 hover:text-white transition-colors"
              disabled={uploading}
            >
              Cancel
            </button>
            <button
              type="submit"
              className={`px-4 py-2 bg-green-600 text-white rounded hover:bg-green-700 transition-colors ${
                uploading ? 'opacity-50 cursor-not-allowed' : ''
              }`}
              disabled={uploading}
            >
              {uploading ? 'Uploading...' : 'Upload'}
            </button>
          </div>
        </form>
      </div>
    </div>
  );
};

export default UploadCsvModal; 
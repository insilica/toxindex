import React, { useEffect, useState } from 'react';
import Modal from 'react-modal';
import { useEnvironment } from '../../context/EnvironmentContext';
import LoadingSpinner from './LoadingSpinner';
import { FaFileCsv, FaFileAlt, FaFileCode, FaDatabase, FaFileImage, FaFile, FaFolderOpen } from 'react-icons/fa';
import { middleEllipsis } from './filename';

interface FilePickerModalProps {
  open: boolean;
  onClose: () => void;
  environmentId?: string | null;
  environments?: any[];
  onPick?: (fileId: string, fileName?: string) => void;
}

// Helper to get file icon based on extension (copied from Layout.tsx)
function getFileIcon(filename: string) {
  const iconClass = "w-5 h-5";
  const ext = filename.split('.').pop()?.toLowerCase();
  if (!ext) return <FaFile className={`text-gray-500 ${iconClass}`} title="File" />;
  if (ext === 'csv') return <FaFileCsv className={`text-green-400 ${iconClass}`} title="CSV file" />;
  if (ext === 'txt') return <FaFileAlt className={`text-gray-400 ${iconClass}`} title="Text file" />;
  if (ext === 'json') return <FaFileCode className={`text-yellow-400 ${iconClass}`} title="JSON file" />;
  if (ext === 'parquet') return <FaDatabase className={`text-blue-400 ${iconClass}`} title="Parquet file" />;
  if (["png","jpg","jpeg","gif","bmp","webp"].includes(ext)) return <FaFileImage className={`text-purple-400 ${iconClass}`} title="Image file" />;
  return <FaFile className={`text-gray-500 ${iconClass}`} title="File" />;
}


const FilePickerModal: React.FC<FilePickerModalProps> = ({ open, onClose, environmentId, environments = [], onPick }) => {
  const [files, setFiles] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const { environments: ctxEnvs } = useEnvironment();

  const envList = environments.length > 0 ? environments : ctxEnvs;
  const env = envList.find((e: any) => e.environment_id === environmentId);

  useEffect(() => {
    if (!open || !environmentId) return;
    setLoading(true);
    setError(null);
    fetch(`/api/environments/${environmentId}/files`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setFiles(data.files || []))
      .catch(() => setError('Failed to load files'))
      .finally(() => setLoading(false));
  }, [open, environmentId]);

  return (
    <Modal
      isOpen={open}
      onRequestClose={onClose}
      contentLabel="Pick a file"
      style={{
        content: { maxWidth: 400, margin: 'auto', minHeight: 200, maxHeight: 550, background: '#181f1b', color: '#fff', borderRadius: 16, overflow: 'auto', top: '0%', left: '62%', right: 'auto', display: 'flex', flexDirection: 'column', padding: 20 },
        overlay: { backgroundColor: 'rgba(0,0,0,0.7)', zIndex: 12 }
      }}
    >
      <div className="flex flex-col gap-4">
        <div className="flex justify-between items-center mb-2">
          <h2 className="text-xl font-bold text-white">Files in {env?.title || 'Environment'}</h2>
          <button onClick={onClose} className="text-gray-400 hover:text-white text-lg font-bold">×</button>
        </div>
        {loading ? (
          <div className="flex justify-center items-center min-h-[200px]"><LoadingSpinner size="medium" /></div>
        ) : error ? (
          <div className="text-red-400 text-center min-h-[200px] flex items-center justify-center">{error}</div>
        ) : files.length === 0 ? (
          <div className="text-gray-400 text-center min-h-[200px] flex items-center justify-center">No files available for this environment.</div>
        ) : (
          <div className="rounded-lg overflow-hidden shadow-lg bg-gray-900">
            <div className="flex items-center gap-3 px-6 py-3 bg-green-950/60 border-b border-green-900">
              <FaFolderOpen className="text-yellow-400 text-xl" />
              <span className="font-bold text-green-200 text-lg truncate" title={env?.title || 'Environment'}>
                {env?.title || 'Environment'}
              </span>
            </div>
            <ul className="divide-y divide-gray-800">
              {files.map(file => {
                return (
                  <li key={file.file_id} className="pl-9 pr-6 py-4 hover:bg-green-900/30 transition flex items-center gap-4 cursor-pointer group rounded-none"
                    onClick={() => onPick?.(file.file_id, file.filename)}
                    style={{ minHeight: 56 }}
                  >
                    <span className="flex items-center justify-center w-8 h-8 rounded-full group-hover:bg-green-800/40 transition bg-gray-800 flex-shrink-0" style={{ aspectRatio: '1 / 1', minWidth: 32, minHeight: 32, maxWidth: 32, maxHeight: 32 }}>
                      {getFileIcon(file.filename)}
                    </span>
                    <span className="font-semibold text-green-200 group-hover:text-green-300 transition flex-1 min-w-0" title={file.filename}>
                      {middleEllipsis(file.filename, 28)}
                    </span>
                    <span className="ml-auto text-xs text-gray-400 font-mono group-hover:text-white transition flex-shrink-0" style={{ minWidth: 110, textAlign: 'right' }}>
                      {file.upload_date ? new Date(file.upload_date).toLocaleString(undefined, { year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' }) : ''}
                    </span>
                  </li>
                );
              })}
            </ul>
          </div>
        )}
      </div>
    </Modal>
  );
};

export default FilePickerModal; 
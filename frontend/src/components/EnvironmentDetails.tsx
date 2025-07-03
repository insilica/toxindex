import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import FilePreviewModal from './FilePreviewModal';
import UploadCsvModal from './UploadCsvModal';
import { FaEye, FaDownload, FaTrash, } from 'react-icons/fa';
import EnvironmentSelector from './shared/EnvironmentSelector';
import LoadingSpinner from './shared/LoadingSpinner';
import { useEnvironment } from "../context/EnvironmentContext";

interface Environment {
  environment_id: string;
  title: string;
  description?: string;
  created_at?: string;
}

interface EnvironmentFile {
  file_id: string;
  filename: string;
  upload_date: string;
}

export const EnvironmentDetails: React.FC = () => {
  const navigate = useNavigate();
  const [env, setEnv] = useState<Environment | null>(null);
  const [files, setFiles] = useState<EnvironmentFile[]>([]);
  const [loading, setLoading] = useState(true);
  const [deleting, setDeleting] = useState(false);
  const [deleteError, setDeleteError] = useState<string | null>(null);
  const [fileDeleteError, setFileDeleteError] = useState<string | null>(null);
  const [previewFileId, setPreviewFileId] = useState<string | null>(null);
  const [previewOpen, setPreviewOpen] = useState(false);
  const [showUploadModal, setShowUploadModal] = useState(false);
  const [tasks, setTasks] = useState<any[]>([]);
  const { selectedEnv, setSelectedEnv, environments, refetchEnvironments } = useEnvironment();

  // Load environment data
  useEffect(() => {
    const loadEnvironment = async () => {
      setLoading(true);
      try {
        if (selectedEnv === "__add__" || selectedEnv === "__manage__") {
          setLoading(false);
          return;
        }

        const envResponse = await fetch(`/api/environments/${selectedEnv}`, { 
          credentials: 'include',
          cache: "no-store"
        });

        const envData = await envResponse.json();
        if (envData.environment) {
          setEnv(envData.environment);
        }
        
        // Fetch files
        const filesResponse = await fetch(`/api/environments/${selectedEnv}/files`, { 
          credentials: 'include',
          cache: "no-store"
        });
        const filesData = await filesResponse.json();
        setFiles(filesData.files || []);

        // Fetch tasks
        const tasksResponse = await fetch(`/api/tasks?environment_id=${selectedEnv}`, { 
          credentials: 'include',
          cache: "no-store"
        });
        const tasksData = await tasksResponse.json();
        setTasks((tasksData.active_tasks || []).concat(tasksData.archived_tasks || []));
      } catch (error) {
        console.error('Error loading environment:', error);
      } finally {
        setLoading(false);
      }
    };

    loadEnvironment();
  }, [selectedEnv, setSelectedEnv]);

  const handleDelete = async () => {
    if (!window.confirm('Are you sure you want to delete this environment? This cannot be undone.')) return;
    setDeleting(true);
    setDeleteError(null);
    try {
      const res = await fetch(`/api/environments/${selectedEnv}`, {
        method: 'DELETE',
        credentials: 'include',
      });
      const data = await res.json();
      if (data.success) {
        // Refresh both files and environments lists
        if (refetchEnvironments) {
          await refetchEnvironments();
        }
        navigate('/settings/environments');
      } else {
        setDeleteError('Failed to delete environment.');
      }
    } catch {
      setDeleteError('Failed to delete environment.');
    } finally {
      setDeleting(false);
    }
  };

  const handleFileDelete = async (file_id: string) => {
    if (!selectedEnv || selectedEnv === "__add__" || selectedEnv === "__manage__") return;
    if (!window.confirm('Are you sure you want to delete this file?')) return;
    setFileDeleteError(null);
    try {
      const res = await fetch(`/api/environments/${selectedEnv}/files/${file_id}`, {
        method: 'DELETE',
        credentials: 'include',
      });
      const data = await res.json();
      if (data.success) {
        setFiles(files => files.filter(f => f.file_id !== file_id));
        // Refresh files list after deletion
        if (refetchEnvironments) {
          await refetchEnvironments();
        }
      } else {
        setFileDeleteError('Failed to delete file.');
      }
    } catch {
      setFileDeleteError('Failed to delete file.');
    }
  };

  if (loading) return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="flex-1 flex items-center justify-center">
        <LoadingSpinner size="large" />
      </div>
    </div>
  );
  
  if (!env) return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="flex-1 flex items-center justify-center">
        <div className="text-white text-xl">Environment not found</div>
      </div>
    </div>
  );

  return (
    <div className="w-full h-full min-h-screen pt-24 px-100 pb-12 text-white flex flex-col" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="flex items-center justify-between mb-2">
        <EnvironmentSelector/>
        <button
          className="px-6 font-semibold transition border border-red-500"
          style={{
            fontWeight: 600,
            fontSize: '.8rem',
            paddingTop: 0,
            paddingBottom: 2,
            lineHeight: 1.1,
            background: 'rgba(220,38,38,0.08)',
            color: '#ef4444',
            borderRadius: '9999px',
            width: 'auto',
            minWidth: 0,
            boxShadow: '0 1px 4px 0 rgba(0,0,0,0.04)',
            borderWidth: 1.5,
            borderStyle: 'solid',
            borderColor: '#ef4444',
            cursor: deleting ? 'not-allowed' : 'pointer',
          }}
          onMouseOver={e => (e.currentTarget.style.background = 'rgba(220,38,38,0.18)')}
          onMouseOut={e => (e.currentTarget.style.background = 'rgba(220,38,38,0.08)')}
          onClick={handleDelete}
          disabled={deleting}
        >
          {deleting ? 'Deleting...' : (
            <span style={{display: 'flex', alignItems: 'center'}}>
              <span style={{fontSize: '1.2rem', fontWeight:900, marginRight:10, lineHeight: 1}}>-</span>
              <span style={{fontSize: '.8rem', fontWeight: 600}}>Delete Environment</span>
            </span>
          )}
        </button>
      </div>
      <div className="w-full border-b border-gray-400 mb-4"></div>
      {env.description && <div className="mb-4 text-gray-300">{env.description}</div>}
      <div className="mb-8 text-sm text-gray-400">
        Created at: {env.created_at && new Date(env.created_at).toLocaleDateString(undefined, { month: 'long', day: 'numeric', year: 'numeric' })}
      </div>
      {deleteError && <div className="text-red-400 mb-4">{deleteError}</div>}

      <div className="pl-8">
        <div className="flex items-center justify-between mb-2 gap-4">
          <h2 className="text-lg font-semibold">Files</h2>
          <div className="flex items-center">
            <button
              className="px-6 font-semibold transition border border-green-600"
              style={{
                fontWeight: 600,
                fontSize: '.8rem',
                paddingTop: 0,
                paddingBottom: 2,
                lineHeight: 1.1,
                background: 'rgba(22,163,74,0.08)',
                color: '#22c55e',
                borderRadius: '9999px',
                width: 'auto',
                minWidth: 0,
                boxShadow: '0 1px 4px 0 rgba(0,0,0,0.04)',
                borderWidth: 1.5,
                borderStyle: 'solid',
                borderColor: '#22c55e',
                cursor: 'pointer',
              }}
              onMouseOver={e => (e.currentTarget.style.background = 'rgba(22,163,74,0.18)')}
              onMouseOut={e => (e.currentTarget.style.background = 'rgba(22,163,74,0.08)')}
              onClick={() => setShowUploadModal(true)}
            >
              <span style={{display: 'flex', alignItems: 'center'}}>
                <span style={{fontSize: '1.2rem', fontWeight:900, marginRight:10, lineHeight: 1}}>+</span>
                <span style={{fontSize: '.8rem', fontWeight: 600}}>Upload CSV</span>
              </span>
            </button>
          </div>
        </div>
        <div className="border-b border-gray-700 mb-4 w-full"></div>
      </div>
      {loading ? (
        <div className="flex justify-center">
          <LoadingSpinner size="medium" text="Loading files..." />
        </div>
      ) : files.length === 0 ? (
        <div className="text-gray-400 mb-6 pl-8">No files uploaded for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700 mb-6 pl-8">
          {files.map(file => (
            <li key={file.file_id} className="py-2 flex items-center justify-between gap-4">
              <div className="flex items-center gap-2">
                <button
                  className="text-left font-medium text-gray-200 hover:text-purple-400 transition truncate max-w-[220px]"
                  style={{ background: 'none', border: 'none', padding: 0, margin: 0, cursor: 'pointer', fontSize: '1rem' }}
                  onClick={() => { setPreviewFileId(file.file_id); setPreviewOpen(true); }}
                  title="Preview file"
                >
                  {file.filename}
                </button>
              </div>
              <div className="flex gap-4 px-3 py-1 rounded-full shadow-sm"
                style={{ background: 'rgba(139, 81, 196, 0.18)', minWidth: 140, justifyContent: 'flex-end' }}>
                <button
                  className="bg-purple-700 hover:bg-purple-600 active:bg-purple-800 text-white transition flex items-center justify-center"
                  style={{ width: 36, height: 36, borderRadius: '50%', fontSize: '1.1rem', display: 'flex', alignItems: 'center', justifyContent: 'center', padding: 0, border: 'none' }}
                  onClick={() => { setPreviewFileId(file.file_id); setPreviewOpen(true); }}
                  title="Preview"
                >
                  <FaEye />
                </button>
                <a
                  href={`/api/environments/${selectedEnv || undefined}/files/${file.file_id}/download`}
                  className="bg-purple-700 hover:bg-purple-600 active:bg-purple-800 text-white transition flex items-center justify-center"
                  style={{ width: 36, height: 36, borderRadius: '50%', fontSize: '1.1rem', display: 'flex', alignItems: 'center', justifyContent: 'center', padding: 0, border: 'none' }}
                  target="_blank"
                  rel="noopener noreferrer"
                  title="Download"
                >
                  <FaDownload />
                </a>
                <button
                  className="bg-purple-700 hover:bg-purple-600 active:bg-purple-800 text-white transition flex items-center justify-center"
                  style={{ width: 36, height: 36, borderRadius: '50%', fontSize: '1.1rem', display: 'flex', alignItems: 'center', justifyContent: 'center', padding: 0, border: 'none' }}
                  onClick={() => handleFileDelete(file.file_id)}
                  title="Delete"
                >
                  <FaTrash />
                </button>
              </div>
            </li>
          ))}
        </ul>
      )}
      {fileDeleteError && <div className="text-red-400 mb-4">{fileDeleteError}</div>}

      <div className="pl-8">
        <h2 className="text-lg font-semibold mb-2">Tasks</h2>
        <div className="border-b border-gray-700 mb-4 w-full"></div>
      </div>
      {tasks.length === 0 ? (
        <div className="text-gray-400 pl-8">No tasks for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700 pl-8">
          {tasks.map(task => (
            <li
              key={task.task_id}
              className="py-2 flex items-center justify-between hover:bg-purple-900/20 rounded-lg cursor-pointer transition-colors"
              style={{ minHeight: 44 }}
              onClick={() => navigate(`/task/${task.task_id}`)}
            >
              <span className="text-white text-left font-medium truncate" style={{ fontSize: '1rem', maxWidth: 220 }} title={task.title}>
                {task.title}
              </span>
              <span
                className="text-xs text-gray-400 ml-2 whitespace-nowrap"
                style={{ fontWeight: 500 }}
                onClick={e => e.stopPropagation()}
              >
                {task.created_at && new Date(task.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric', year: 'numeric' })}
              </span>
            </li>
          ))}
        </ul>
      )}

      {(selectedEnv) && (
        <FilePreviewModal
          fileId={previewFileId}
          isOpen={previewOpen}
          onRequestClose={() => setPreviewOpen(false)}
        />
      )}
      
      {showUploadModal && (
        <UploadCsvModal
          open={showUploadModal}
          onClose={() => setShowUploadModal(false)}
          environments={environments}
          defaultEnvId={selectedEnv || undefined}
          onUploadSuccess={() => {
            setShowUploadModal(false);
            if (selectedEnv && selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
              fetch(`/api/environments/${selectedEnv}/files`, { credentials: "include", cache: "no-store" })
                .then(res => res.json())
                .then(data => setFiles(data.files || []));
            }
          }}
        />
      )}
    </div>
  );
};

export default EnvironmentDetails;
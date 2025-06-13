import React, { useState, useEffect } from "react";
import { useNavigate, useParams } from "react-router-dom";
import FilePreviewModal from './FilePreviewModal';
import { FaEye, FaDownload, FaTrash, FaFileCsv, FaFileAlt, FaFileCode, FaDatabase, FaFileImage, FaFile } from 'react-icons/fa';
import UploadCsvModal from './UploadCsvModal';

const CreateEnvironment: React.FC = () => {
  const [title, setTitle] = useState("");
  const [description, setDescription] = useState("");
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);

    const res = await fetch("/api/environment/new", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      credentials: "include",
      body: JSON.stringify({ title, description }),
    });

    if (res.ok) {
      const data = await res.json();
      if (data.environment_id) {
        navigate(`/environment/${data.environment_id}`);
      }
    } else {
      setError("Failed to create environment.");
    }
  };

  return (
    <form onSubmit={handleSubmit} className="max-w-md mx-auto p-6 bg-white rounded shadow">
      <h2 className="text-xl font-bold mb-4">Create Environment</h2>
      <input
        type="text"
        placeholder="Environment name"
        value={title}
        onChange={e => setTitle(e.target.value)}
        required
        className="w-full mb-3 p-2 border rounded"
      />
      <textarea
        placeholder="Description (optional)"
        value={description}
        onChange={e => setDescription(e.target.value)}
        className="w-full mb-3 p-2 border rounded"
      />
      {error && <div className="text-red-500 mb-2">{error}</div>}
      <button
        type="submit"
        className="w-full bg-blue-600 text-white py-2 rounded hover:bg-blue-700"
      >
        Create
      </button>
    </form>
  );
};

interface FileMeta {
  file_id: number;
  filename: string;
  filepath: string;
  created_at: string;
}

interface EnvDetails {
  environment_id: string;
  title: string;
  description?: string;
  created_at?: string;
}

export const EnvironmentDetails: React.FC<{ refreshEnvFiles?: () => void }> = ({ refreshEnvFiles }) => {
  const { env_id } = useParams<{ env_id: string }>();
  const [env, setEnv] = useState<EnvDetails | null>(null);
  const [environments, setEnvironments] = useState<EnvDetails[]>([]);
  const [files, setFiles] = useState<FileMeta[]>([]);
  const [loading, setLoading] = useState(true);
  const [deleting, setDeleting] = useState(false);
  const [deleteError, setDeleteError] = useState<string | null>(null);
  const [tasks, setTasks] = useState<any[]>([]);
  const [fileDeleteError, setFileDeleteError] = useState<string | null>(null);
  const [uploading, setUploading] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const fileInputRef = React.useRef<HTMLInputElement>(null);
  const navigate = useNavigate();
  const [previewFileId, setPreviewFileId] = useState<number | null>(null);
  const [previewOpen, setPreviewOpen] = useState(false);
  const [showUploadModal, setShowUploadModal] = useState(false);

  useEffect(() => {
    if (!env_id) return;
    setLoading(true);
    fetch(`/api/environments`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setEnvironments(data.environments || []);
        const found = (data.environments || []).find((e: any) => e.environment_id == env_id);
        setEnv(found || null);
      });
    fetch(`/api/environment/${env_id}/files`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setFiles(data.files || []))
      .finally(() => setLoading(false));
  }, [env_id]);

  useEffect(() => {
    if (!env_id) return;
    fetch(`/tasks?environment_id=${env_id}`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setTasks((data.active_tasks || []).concat(data.archived_tasks || [])));
  }, [env_id]);

  const handleDelete = async () => {
    if (!window.confirm('Are you sure you want to delete this environment? This cannot be undone.')) return;
    setDeleting(true);
    setDeleteError(null);
    try {
      const res = await fetch(`/environment/${env_id}`, {
        method: 'DELETE',
        credentials: 'include',
      });
      const data = await res.json();
      if (data.success) {
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

  const handleFileDelete = async (file_id: number) => {
    if (!window.confirm('Are you sure you want to delete this file?')) return;
    setFileDeleteError(null);
    try {
      const res = await fetch(`/api/file/${file_id}`, {
        method: 'DELETE',
        credentials: 'include',
      });
      const data = await res.json();
      if (data.success) {
        setFiles(files => files.filter(f => f.file_id !== file_id));
      } else {
        setFileDeleteError('Failed to delete file.');
      }
    } catch {
      setFileDeleteError('Failed to delete file.');
    }
  };

  const handleUploadClick = () => {
    if (fileInputRef.current) fileInputRef.current.value = '';
    fileInputRef.current?.click();
  };

  const handleFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setUploading(true);
    setUploadError(null);
    const formData = new FormData();
    formData.append('file', file);
    formData.append('environment_id', env_id!);
    try {
      const res = await fetch('/api/upload-file', {
        method: 'POST',
        credentials: 'include',
        body: formData,
      });
      if (!res.ok) throw new Error('Upload failed');
      // Refresh file list
      fetch(`/api/environment/${env_id}/files`, { credentials: 'include' })
        .then(res => res.json())
        .then(data => setFiles(data.files || []));
      if (refreshEnvFiles) refreshEnvFiles();
    } catch {
      setUploadError('Failed to upload file.');
    } finally {
      setUploading(false);
    }
  };

  function getFileIcon(filename: string) {
    const ext = filename.split('.').pop()?.toLowerCase();
    if (!ext) return <FaFile />;
    if (ext === 'csv') return <FaFileCsv className="text-green-400" title="CSV file" />;
    if (ext === 'txt') return <FaFileAlt className="text-gray-400" title="Text file" />;
    if (ext === 'json') return <FaFileCode className="text-yellow-400" title="JSON file" />;
    if (ext === 'parquet') return <FaDatabase className="text-blue-400" title="Parquet file" />;
    if (["png","jpg","jpeg","gif","bmp","webp"].includes(ext)) return <FaFileImage className="text-purple-400" title="Image file" />;
    return <FaFile className="text-gray-500" title="File" />;
  }

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!env) return <div className="text-white p-8">Environment not found.</div>;

  return (
    <div className="w-full h-full min-h-screen pt-24 px-100 pb-12 text-white flex flex-col" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="flex items-center justify-between mb-2">
        <div className="relative group" style={{ minWidth: 180, maxWidth: 340 }}>
          <select
            className="font-bold text-white text-lg px-0 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-gray-900 group-hover:border group-hover:border-gray-700 group-hover:px-6 group-hover:pr-10 group-hover:cursor-pointer focus:bg-gray-900 focus:border focus:border-gray-700 focus:px-6 focus:pr-10"
            style={{ fontSize: '1.2rem', minWidth: 180, maxWidth: 340, cursor: 'pointer' }}
            value={env_id}
            onChange={e => navigate(`/environment/${e.target.value}`)}
          >
            {environments.map(e => (
              <option key={e.environment_id} value={e.environment_id}>{e.title}</option>
            ))}
          </select>
          {/* Chevron icon, only visible on hover/focus */}
          <span
            className="pointer-events-none absolute right-2 top-1/2 transform -translate-y-1/2 text-gray-400 group-hover:opacity-100 group-focus-within:opacity-100 opacity-0 transition-opacity duration-150"
            style={{ right: '1.2rem' }}
          >
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M7 10l5 5 5-5" />
            </svg>
          </span>
        </div>
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
          {deleting ? 'Deleting...' : 
          <span style={{display: 'flex', alignItems: 'center'}}>
            <span style={{fontSize: '1.2rem', fontWeight:900, marginRight:10, lineHeight: 1}}>-</span>
            <span style={{fontSize: '.8rem', fontWeight: 600}}>Delete Environment</span>
          </span>}
        </button>
      </div>
      <div className="w-full border-b border-gray-400 mb-4"></div>
      {env.description && <div className="mb-4 text-gray-300">{env.description}</div>}
      <div className="mb-40 text-sm text-gray-400">
        Created at: {env.created_at && new Date(env.created_at).toLocaleDateString(undefined, { month: 'long', day: 'numeric', year: 'numeric' })}
      </div>
      {deleteError && <div className="text-red-400 mb-4">{deleteError}</div>}
      <div className="pl-8">
        <div className="flex items-center justify-between mb-2 gap-4">
          <h2 className="text-lg font-semibold">Uploaded Files</h2>
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
                cursor: uploading ? 'not-allowed' : 'pointer',
              }}
              onMouseOver={e => (e.currentTarget.style.background = 'rgba(22,163,74,0.18)')}
              onMouseOut={e => (e.currentTarget.style.background = 'rgba(22,163,74,0.08)')}
              onClick={() => setShowUploadModal(true)}
              disabled={uploading}
            >
              {uploading ? 'Uploading...' : (
                <span style={{display: 'flex', alignItems: 'center'}}>
                  <span style={{fontSize: '1.2rem', fontWeight:900, marginRight:10, lineHeight: 1}}>+</span>
                  <span style={{fontSize: '.8rem', fontWeight: 600}}>Upload CSV</span>
                </span>
              )}
            </button>
          </div>
        </div>
        <div className="border-b border-gray-700 mb-4 w-full"></div>
      </div>
      {files.length === 0 ? (
        <div className="text-gray-400 mb-6 pl-8">No files uploaded for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700 mb-6 pl-8">
          {files.map(file => (
            <li key={file.file_id} className="py-2 flex items-center justify-between gap-4">
              <div className="flex items-center gap-2">
                {getFileIcon(file.filename)}
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
                  href={`/api/file/${file.file_id}/download`}
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
      {uploadError && <div className="text-red-400 mb-2">{uploadError}</div>}
      <div className="pl-8">
        <h2 className="text-lg font-semibold mb-2 mt-8">List of Tasks</h2>
        <div className="border-b border-gray-700 mb-4 w-full"></div>
      </div>
      {tasks.length === 0 ? (
        <div className="text-gray-400 pl-8">No tasks for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700 pl-8">
          {tasks.map(task => (
            <li key={task.task_id} className="py-2 flex items-center justify-between">
              <span>{task.title}</span>
              <span className="text-xs text-gray-400 ml-2">{task.created_at && new Date(task.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric', year: 'numeric' })}</span>
              {task.archived && <span className="ml-2 text-yellow-400 text-xs">(archived)</span>}
            </li>
          ))}
        </ul>
      )}
      <FilePreviewModal
        fileId={previewFileId}
        isOpen={previewOpen}
        onRequestClose={() => setPreviewOpen(false)}
      />
      <UploadCsvModal
        open={showUploadModal}
        onClose={() => setShowUploadModal(false)}
        environments={environments}
        defaultEnvId={env_id}
        onUploadSuccess={() => {
          setShowUploadModal(false);
          // Refresh file list after upload
          fetch(`/api/environment/${env_id}/files`, { credentials: 'include' })
            .then(res => res.json())
            .then(data => setFiles(data.files || []));
          if (refreshEnvFiles) refreshEnvFiles();
        }}
      />
    </div>
  );
};

export default CreateEnvironment; 
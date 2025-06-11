import React, { useState, useEffect } from "react";
import { useNavigate, useParams } from "react-router-dom";

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

export const EnvironmentDetails: React.FC = () => {
  const { env_id } = useParams<{ env_id: string }>();
  const [env, setEnv] = useState<EnvDetails | null>(null);
  const [files, setFiles] = useState<FileMeta[]>([]);
  const [loading, setLoading] = useState(true);
  const [deleting, setDeleting] = useState(false);
  const [deleteError, setDeleteError] = useState<string | null>(null);
  const [tasks, setTasks] = useState<any[]>([]);
  const [fileDeleteError, setFileDeleteError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    if (!env_id) return;
    setLoading(true);
    fetch(`/api/environments`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
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

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!env) return <div className="text-white p-8">Environment not found.</div>;

  return (
    <div className="w-full h-full min-h-screen pt-24 px-100 pb-12 text-white flex flex-col" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="flex items-center justify-between mb-2">
        <h1 className="font-bold text-white" style={{ fontSize: '1rem' }}>{env.title}</h1>
        <button
          className="px-6 font-semibold transition"
          style={{
            fontWeight: 600,
            fontSize: '.8rem',
            paddingTop: 0,
            paddingBottom: 0,
            lineHeight: 1.1,
            boxShadow: '0 1px 4px 0 rgba(0,0,0,0.04)',
            background: '#f3f4f6',
            color: '#111',
            borderRadius: '9999px',
            width: 'auto',
            minWidth: 0,
          }}
          onMouseOver={e => (e.currentTarget.style.background = '#e5e7eb')}
          onMouseOut={e => (e.currentTarget.style.background = '#f3f4f6')}
          onClick={handleDelete}
          disabled={deleting}
        >
          {deleting ? 'Deleting...' : 'Delete Environment'}
        </button>
      </div>
      <div className="w-full border-b border-gray-400 mb-4"></div>
      {env.description && <div className="mb-4 text-gray-300">{env.description}</div>}
      <div className="mb-6 text-sm text-gray-400">
        Created at: {env.created_at && new Date(env.created_at).toLocaleDateString(undefined, { month: 'long', day: 'numeric', year: 'numeric' })}
      </div>
      {deleteError && <div className="text-red-400 mb-4">{deleteError}</div>}
      <h2 className="text-lg font-semibold mb-2">Uploaded Files</h2>
      <div className="w-full border-b border-gray-700 mb-4"></div>
      {files.length === 0 ? (
        <div className="text-gray-400 mb-6">No files uploaded for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700 mb-6">
          {files.map(file => (
            <li key={file.file_id} className="py-2 flex items-center justify-between gap-4">
              <span>{file.filename}</span>
              <div className="flex gap-2">
                <a
                  href={`/api/file/${file.file_id}/download`}
                  className="text-blue-400 hover:underline px-3 py-1 rounded bg-gray-800 hover:bg-gray-700"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  Download
                </a>
                <button
                  className="px-3 py-1 rounded bg-red-700 hover:bg-red-800 text-white text-sm"
                  onClick={() => handleFileDelete(file.file_id)}
                >
                  Delete
                </button>
              </div>
            </li>
          ))}
        </ul>
      )}
      {fileDeleteError && <div className="text-red-400 mb-4">{fileDeleteError}</div>}
      <h2 className="text-lg font-semibold mb-2 mt-8">List of Tasks</h2>
      <div className="w-full border-b border-gray-700 mb-4"></div>
      {tasks.length === 0 ? (
        <div className="text-gray-400">No tasks for this environment.</div>
      ) : (
        <ul className="divide-y divide-gray-700">
          {tasks.map(task => (
            <li key={task.task_id} className="py-2 flex items-center justify-between">
              <span>{task.title}</span>
              <span className="text-xs text-gray-400 ml-2">{task.created_at && new Date(task.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric', year: 'numeric' })}</span>
              {task.archived && <span className="ml-2 text-yellow-400 text-xs">(archived)</span>}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
};

export default CreateEnvironment; 
import React, { useEffect, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import FilePreviewModal from './FilePreviewModal';

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
  const [, setEnvironments] = useState<EnvDetails[]>([]);
  const [, setFiles] = useState<FileMeta[]>([]);
  const [loading, setLoading] = useState(true);
  const [, setDeleting] = useState(false);
  const [, setDeleteError] = useState<string | null>(null);
  const navigate = useNavigate();
  const [previewFileId, ] = useState<number | null>(null);
  const [previewOpen, setPreviewOpen] = useState(false);
  const [, setTasks] = useState<any[]>([]);
  // const [, setFileDeleteError] = useState<string | null>(null);

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
    fetch(`/api/environments/${env_id}/files`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setFiles(data.files || []))
      .finally(() => setLoading(false));
  }, [env_id]);

  useEffect(() => {
    if (!env_id) return;
    fetch(`/api/tasks?environment_id=${env_id}`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setTasks((data.active_tasks || []).concat(data.archived_tasks || [])));
  }, [env_id]);

  const handleDelete = async () => {
    if (!window.confirm('Are you sure you want to delete this environment? This cannot be undone.')) return;
    setDeleting(true);
    setDeleteError(null);
    try {
      const res = await fetch(`/api/environments/${env_id}`, {
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

  // const handleFileDelete = async (file_id: number) => {
  //   if (!window.confirm('Are you sure you want to delete this file?')) return;
  //   try {
  //     const res = await fetch(`/api/environments/${env_id}/files/${file_id}`, {
  //       method: 'DELETE',
  //       credentials: 'include',
  //     });
  //     const data = await res.json();
  //     if (data.success) {
  //       setFiles(files => files.filter(f => f.file_id !== file_id));
  //     } else {
  //       setFileDeleteError('Failed to delete file.');
  //     }
  //   } catch {
  //     setFileDeleteError('Failed to delete file.');
  //   }
  // };

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!env) return <div className="text-white p-8">Environment not found.</div>;

  return (
    <div className="w-full h-full min-h-screen pt-24 px-100 pb-12 text-white flex flex-col" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <FilePreviewModal
        fileId={previewFileId}
        envId={env_id || ''}
        isOpen={previewOpen}
        onRequestClose={() => setPreviewOpen(false)}
      />
    </div>
  );
};

export default EnvironmentDetails;
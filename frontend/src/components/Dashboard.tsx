import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { FaListAlt, FaPlus } from 'react-icons/fa';

interface Environment {
  environment_id: string;
  title: string;
}

interface Task {
  task_id: string;
  title: string;
  archived: boolean;
  // Add other fields as needed
}

const Dashboard: React.FC = () => {
  const [environments, setEnvironments] = useState<Environment[]>([]);
  const [selectedEnv, setSelectedEnv] = useState<string>("");
  const [chatInput, setChatInput] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [user, setUser] = useState<{ email?: string } | null>(null);
  const [tasksLoading, setTasksLoading] = useState(false);
  const [activeTasks, setActiveTasks] = useState<Task[]>([]);
  const [archivedTasks, setArchivedTasks] = useState<Task[]>([]);
  const [activeTab, setActiveTab] = useState<'tasks' | 'archive'>('tasks');
  const [uploading, setUploading] = useState(false);
  const fileInputRef = React.useRef<HTMLInputElement>(null);
  const navigate = useNavigate();
  const [showUploadModal, setShowUploadModal] = useState(false);
  const [uploadEnvId, setUploadEnvId] = useState<string>("");
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [uploadError, setUploadError] = useState<string | null>(null);

  useEffect(() => {
    fetch("/environments", { credentials: "include" })
      .then(res => res.json())
      .then(data => {
        console.log("Fetched environments:", data.environments);
        setEnvironments(data.environments || []);
        if (data.environments && data.environments.length > 0) {
          setSelectedEnv(data.environments[0].environment_id);
        } else {
          setSelectedEnv("__add__");
        }
      });
  }, []);

  useEffect(() => {
    fetch("/api/me", { credentials: "include", cache: "no-store" })
      .then(res => res.ok ? res.json() : null)
      .then(data => setUser(data));
  }, []);

  useEffect(() => {
    setTasksLoading(true);
    fetch("/tasks", { credentials: "include" })
      .then(res => res.json())
      .then(data => {
        setActiveTasks(data.active_tasks || []);
        setArchivedTasks(data.archived_tasks || []);
      })
      .finally(() => setTasksLoading(false));
  }, []);

  const handleChatSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    // You can handle chat input here (send to backend, etc.)
    setChatInput("");
  };

  const handleDropdownChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    if (e.target.value === "__manage__") {
      navigate("/settings/environments");
    } else if (e.target.value === "__add__") {
      navigate("/settings/environments/create");
    } else {
      setSelectedEnv(e.target.value);
    }
  };

  const handleLogout = () => {
    fetch("/api/logout", {
      method: "GET",
      credentials: "include",
    }).then(() => {
      window.location.href = "/login";
    });
  };

  const archiveTask = (task_id: string) => {
    fetch(`/task/${task_id}/archive`, { method: "POST", credentials: "include" })
      .then(res => res.json())
      .then(() => {
        setActiveTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const archived = activeTasks.find(t => t.task_id === task_id);
        if (archived) setArchivedTasks(tasks => [archived, ...tasks]);
      });
  };

  const unarchiveTask = (task_id: string) => {
    fetch(`/task/${task_id}/unarchive`, { method: "POST", credentials: "include" })
      .then(res => res.json())
      .then(() => {
        setArchivedTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const unarchived = archivedTasks.find(t => t.task_id === task_id);
        if (unarchived) setActiveTasks(tasks => [unarchived, ...tasks]);
      });
  };

  const handlePlusClick = () => {
    setUploadEnvId(selectedEnv);
    setShowUploadModal(true);
    setUploadFile(null);
    setUploadError(null);
  };

  const handleUploadConfirm = async () => {
    if (!uploadFile) {
      setUploadError("Please select a file.");
      return;
    }
    setUploading(true);
    setUploadError(null);
    const formData = new FormData();
    formData.append('file', uploadFile);
    formData.append('environment_id', uploadEnvId);
    try {
      const res = await fetch('/api/upload-file', {
        method: 'POST',
        credentials: 'include',
        body: formData,
      });
      if (!res.ok) throw new Error('Upload failed');
      setShowUploadModal(false);
    } catch (err) {
      setUploadError('Failed to upload file.');
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)', position: 'relative' }}>
      <div className="w-full flex justify-center pt-24 pb-4">
        <h2 className="text-2xl font-normal text-white text-center" style={{ fontFamily: 'inherit' }}>Is it toxic, or just misunderstood? Let's break it down!</h2>
      </div>
      {/* Tabs and Task Lists */}
      <div className="w-full max-w-2xl mx-auto mt-8">
        <div className="flex space-x-2 mb-4 border-b border-gray-700">
          <button
            className={`px-4 py-2 font-semibold flex items-center gap-2 focus:outline-none transition border-b-4 ${activeTab === 'tasks' ? 'border-white text-white' : 'border-transparent text-gray-400 hover:text-white'}`}
            style={{ background: 'none', borderRadius: 0 }}
            onClick={() => setActiveTab('tasks')}
          >
            <FaListAlt className="inline-block" /> Tasks
          </button>
          <button
            className={`px-4 py-2 font-semibold flex items-center gap-2 focus:outline-none transition border-b-4 ${activeTab === 'archive' ? 'border-white text-white' : 'border-transparent text-gray-400 hover:text-white'}`}
            style={{ background: 'none', borderRadius: 0 }}
            onClick={() => setActiveTab('archive')}
          >
            <FaListAlt className="inline-block" /> Archive
          </button>
        </div>
        {activeTab === 'tasks' && (
          tasksLoading ? (
            <div className="text-white">Loading tasks...</div>
          ) : activeTasks.length === 0 ? (
            <div className="text-gray-400">No active tasks.</div>
          ) : (
            <ul className="bg-gray-900 bg-opacity-60 rounded-b-lg shadow divide-y divide-gray-800">
              {activeTasks.map(task => (
                <li key={task.task_id} className="flex items-center justify-between px-4 py-3">
                  <span className="text-white">{task.title}</span>
                  <button
                    onClick={() => archiveTask(task.task_id)}
                    className="ml-4 px-3 py-1 rounded bg-yellow-600 text-white hover:bg-yellow-700 text-sm"
                    title="Archive task"
                  >
                    Archive
                  </button>
                </li>
              ))}
            </ul>
          )
        )}
        {activeTab === 'archive' && (
          tasksLoading ? (
            <div className="text-white">Loading tasks...</div>
          ) : archivedTasks.length === 0 ? (
            <div className="text-gray-400">No archived tasks.</div>
          ) : (
            <ul className="bg-gray-900 bg-opacity-40 rounded-b-lg shadow divide-y divide-gray-800">
              {archivedTasks.map(task => (
                <li key={task.task_id} className="flex items-center justify-between px-4 py-3">
                  <span className="text-gray-300">{task.title}</span>
                  <button
                    onClick={() => unarchiveTask(task.task_id)}
                    className="ml-4 px-3 py-1 rounded bg-green-700 text-white hover:bg-green-800 text-sm"
                    title="Unarchive task"
                  >
                    Unarchive
                  </button>
                </li>
              ))}
            </ul>
          )
        )}
      </div>
      {/* Chatbox and rest of dashboard content below */}
      <div className="flex-1 flex flex-col justify-end">
        <div className="w-full flex justify-center pb-8">
          <form
            onSubmit={handleChatSubmit}
            className="flex flex-col items-center w-full"
            style={{ maxWidth: '900px' }}
          >
            <div className="relative w-full" style={{ maxWidth: '900px', width: '100%' }}>
              <textarea
                rows={3}
                placeholder="Ask me your toxicology question [ Is Aspirin hepatotoxic? ]"
                value={chatInput}
                onChange={e => setChatInput(e.target.value)}
                className="w-full pt-4 pb-4 pr-28 pl-8 text-lg rounded-2xl border border-gray-700 bg-gray-900 bg-opacity-70 text-white resize-none min-h-[80px] shadow-2xl focus:outline-none focus:ring-2 focus:ring-green-400 text-left placeholder:text-left"
                style={{ minHeight: 80, fontFamily: 'inherit', width: '100%', boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)' }}
              />
              <div className="absolute right-4 bottom-11 flex items-center space-x-2 z-30">
                <button
                  className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                  style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                  title="Upload CSV file"
                  onClick={handlePlusClick}
                  disabled={uploading}
                >
                  <FaPlus className="w-5 h-5 text-black" />
                </button>
                <button
                  type="submit"
                  className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                  style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                  title="Submit"
                >
                  <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={2.2} stroke="black" className="w-5 h-5">
                    <path strokeLinecap="round" strokeLinejoin="round" d="M5 10l7-7m0 0l7 7m-7-7v18" />
                  </svg>
                </button>

              </div>
              <div className="absolute left-4 z-20 flex items-end" style={{ position: 'relative', width: '210px', overflow: 'visible', bottom: '3rem' }}>
                <div style={{ position: 'relative', width: '100%' }}>
                  <select
                    className="pl-4 pr-10 bg-black bg-opacity-60 text-white text-sm border border-gray-700 shadow-sm appearance-none focus:ring-2 focus:ring-green-400 w-full"
                    style={{ height: '1.7rem', minHeight: '1.7rem', borderRadius: '999px', lineHeight: '1.2rem', position: 'relative', zIndex: 1, minWidth: '180px', maxWidth: '210px', backgroundColor: 'rgba(16,16,16,0.6)' }}
                    value={selectedEnv}
                    onChange={handleDropdownChange}
                  >
                    {environments.length > 0 && environments.map(env => (
                      <option key={env.environment_id} value={env.environment_id}>
                        {env.title}
                      </option>
                    ))}
                    <option value="__add__">+ Add environment</option>
                    <option value="__manage__">&#9881; Manage environments</option>
                  </select>
                  <span style={{ pointerEvents: 'none', position: 'absolute', right: '10px', top: '50%', transform: 'translateY(-50%)', display: 'flex', alignItems: 'center', height: '100%', zIndex: 2, background: 'transparent', paddingRight: '2px' }}>
                    <svg width="15" height="15" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
                      <path d="M7 10l5 5 5-5" stroke="#fff" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round"/>
                    </svg>
                  </span>
                </div>
              </div>
            </div>
            {error && <div className="text-red-500 mt-2">{error}</div>}
          </form>
        </div>
      </div>
      {/* Upload Modal */}
      {showUploadModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black bg-opacity-60">
          <div className="bg-gray-900 rounded-lg p-6 w-full max-w-md relative">
            <button
              className="absolute top-2 right-2 text-gray-400 hover:text-white"
              onClick={() => setShowUploadModal(false)}
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
              {environments.map(env => (
                <option key={env.environment_id} value={env.environment_id}>{env.title}</option>
              ))}
            </select>
            <label className="block text-white mb-2">CSV File</label>
            <input
              type="file"
              accept=".csv"
              className="w-full mb-4"
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
      )}
    </div>
  );
};

export default Dashboard; 
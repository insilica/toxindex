import React, { useEffect, useState, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { FaListAlt, FaPlus, FaArchive } from 'react-icons/fa';

interface Environment {
  environment_id: string;
  title: string;
}

interface Task {
  task_id: string;
  title: string;
  archived: boolean;
  environment_id: string;
  created_at?: string;
  session_id?: string;
  // Add other fields as needed
}

interface DashboardProps {
  selectedModel?: string;
  selectedEnv?: string;
  setSelectedEnv?: (envId: string) => void;
  environments: Environment[];
  refetchChatSessions?: () => void;
  refetchEnvironments: () => void;
  loadingEnvironments: boolean;
}

const Dashboard: React.FC<DashboardProps> = ({ selectedModel, selectedEnv, setSelectedEnv, environments, refetchChatSessions, loadingEnvironments }) => {
  const [chatInput, setChatInput] = useState("");
  const [error, setError] = useState<string | null>(null);
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
  const [dragActive, setDragActive] = useState(false);
  const [selectWidth, setSelectWidth] = useState(120); // default min width
  const selectRef = useRef<HTMLSelectElement>(null);
  const spanRef = useRef<HTMLSpanElement>(null);

  useEffect(() => {
    setTasksLoading(true);
    let url = "/api/tasks";
    if (selectedEnv && selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
      url += `?environment_id=${selectedEnv}`;
    }
    fetch(url, { credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(data => {
        setActiveTasks(data.active_tasks || []);
        setArchivedTasks(data.archived_tasks || []);
      })
      .finally(() => setTasksLoading(false));
  }, [selectedEnv]);

  useEffect(() => {
    if (spanRef.current) {
      setSelectWidth(spanRef.current.offsetWidth + 40); // add some padding
    }
  }, [selectedEnv, environments]);

  // Auto-select first environment and restore from localStorage
  useEffect(() => {
    if (!setSelectedEnv) return;
    const stored = localStorage.getItem('selectedEnv');
    if (stored && environments.some(e => e.environment_id === stored)) {
      setSelectedEnv(stored);
    } else if (environments.length > 0) {
      setSelectedEnv(environments[0].environment_id);
    }
  }, [environments, setSelectedEnv]);

  useEffect(() => {
    if (selectedEnv) {
      localStorage.setItem('selectedEnv', selectedEnv);
    }
  }, [selectedEnv]);

  const handleDropdownChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    if (e.target.value === "__manage__") {
      navigate("/settings/environments");
    } else if (e.target.value === "__add__") {
      navigate("/settings/environments/create");
    } else {
      setSelectedEnv && setSelectedEnv(String(e.target.value || ''));
    }
  };

  const archiveTask = (task_id: string) => {
    fetch(`/api/tasks/${task_id}/archive`, { method: "POST", credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(() => {
        setActiveTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const archived = activeTasks.find(t => t.task_id === task_id);
        if (archived) setArchivedTasks(tasks => [archived, ...tasks]);
      });
  };

  const unarchiveTask = (task_id: string) => {
    fetch(`/api/tasks/${task_id}/unarchive`, { method: "POST", credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(() => {
        setArchivedTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const unarchived = archivedTasks.find(t => t.task_id === task_id);
        if (unarchived) setActiveTasks(tasks => [unarchived, ...tasks]);
      });
  };

  const handlePlusClick = () => {
    setUploadEnvId(selectedEnv || '');
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
        cache: 'no-store',
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

  const handleFormSubmit = async (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    setError(null);
    if (!chatInput.trim()) {
      setError("Please enter a question.");
      return;
    }
    if (!selectedEnv || selectedEnv === "__add__" || selectedEnv === "__manage__") {
      setError("Please select an environment.");
      return;
    }
    let endpoint = null;
    if (selectedModel === "toxindex-rap") {
      endpoint = "/api/run-probra-task";
    } else if (selectedModel === "toxindex-vanilla") {
      endpoint = "/api/run-vanilla-task";
    } else if (selectedModel === "toxindex-pathway") {
      setError("ToxIndex Pathway is not yet supported.");
      return;
    } else if (selectedModel === "toxindex-4th") {
      setError("ToxIndex 4th is not yet supported.");
      return;
    } else if (selectedModel === "toxindex-5th") {
      setError("ToxIndex 5th is not yet supported.");
      return;
    } else {
      setError("Unknown model selected.");
      return;
    }
    try {
      const res = await fetch(endpoint, {
        method: "POST",
        credentials: "include",
        cache: "no-store",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          prompt: chatInput,
          environment_id: selectedEnv,
        }),
      });
      if (!res.ok) throw new Error("Failed to create task");
      const data = await res.json();
      setChatInput("");
      // Optionally, you can refresh the tasks list here by refetching
      setTasksLoading(true);
      let url = "/api/tasks";
      if (selectedEnv && selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
        url += `?environment_id=${selectedEnv}`;
      }
      const tasksData = await fetch(url, { credentials: "include", cache: "no-store" }).then(r => r.json());
      setActiveTasks(tasksData.active_tasks || []);
      setArchivedTasks(tasksData.archived_tasks || []);
      // --- NEW: Use session_id from backend response for redirect ---
      if (data.session_id) {
        navigate(`/chat/session/${data.session_id}`);
        if (typeof refetchChatSessions === 'function') {
          refetchChatSessions();
        }
      }
    } catch (err) {
      setError("Failed to create new chat session.");
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)', position: 'relative' }}>
      <div className="w-full flex justify-center pt-24 pb-4">
        <h2 className="text-2xl font-normal text-white text-center" style={{ fontFamily: 'inherit' }}>Is it toxic, or just misunderstood? Let's break it down!</h2>
      </div>
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
            className={`px-4 py-0 font-semibold flex items-center gap-2 focus:outline-none transition border-b-4 ${activeTab === 'archive' ? 'border-white text-white' : 'border-transparent text-gray-400 hover:text-white'}`}
            style={{ background: 'none', borderRadius: 0 }}
            onClick={() => setActiveTab('archive')}
          >
            <FaArchive className="inline-block" /> Archive
          </button>
        </div>
        {activeTab === 'tasks' && (
          tasksLoading ? (
            <div className="text-white">Loading tasks...</div>
          ) : activeTasks.length === 0 ? (
            <div className="text-gray-400">No active tasks.</div>
          ) : (
            <ul className="divide-y divide-gray-800">
              {activeTasks.map(task => (
                <li
                  key={task.task_id}
                  className="flex items-center justify-between px-4 py-1 transition-colors duration-150 hover:bg-purple-900/20 rounded-lg text-base"
                  style={{ background: 'none', minHeight: 44 }}
                  onClick={() => navigate(`/task/${task.task_id}`)}
                >
                  <div className="flex items-center min-w-0 flex-1 gap-3">
                    <span
                      className="text-white text-left font-medium truncate"
                      style={{
                        background: 'none',
                        fontSize: '.95rem',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        whiteSpace: 'nowrap',
                        lineHeight: 1.1,
                        boxShadow: 'none',
                        outline: 'none',
                        border: 'none',
                        maxWidth: 180
                      }}
                      title={task.title}
                    >
                      {task.title}
                    </span>
                    {task.created_at && (
                      <span className="text-xs text-gray-400 ml-2 whitespace-nowrap" style={{fontWeight: 500}}>
                        {new Date(task.created_at).toLocaleString('en-US', {
                          month: 'short', day: '2-digit', year: 'numeric',
                          hour: 'numeric', minute: '2-digit', hour12: true
                        }).replace(',', '')}
                      </span>
                    )}
                    {task.session_id && (
                      <span className="ml-2 px-2 py-0.5 bg-gray-800 rounded text-green-400 font-mono text-xs" title={task.session_id} style={{fontWeight: 500}}>
                        chat {task.session_id.slice(0, 6)}
                      </span>
                    )}
                  </div>
                  <button
                    onClick={(e) => {
                      e.stopPropagation();
                      archiveTask(task.task_id);
                    }}
                    className="ml-2 flex items-center justify-center w-7 h-7 rounded-full bg-purple-800 hover:bg-purple-700 active:bg-purple-900 text-white transition border-none shadow-sm focus:outline-none focus:ring-2 focus:ring-purple-400 hover:ring-2 hover:ring-purple-400"
                    style={{ padding: 0, borderRadius: '50%', fontSize: '1.1rem', border: 'none', outline: 'none', boxShadow: 'none', cursor: 'pointer' }}
                    title="Archive task"
                  >
                    <FaArchive />
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
            <ul className="divide-y divide-gray-800">
              {archivedTasks.map(task => (
                <li
                  key={task.task_id}
                  className="flex items-center justify-between px-4 py-3 hover:bg-gray-800 cursor-pointer rounded"
                  onClick={(e) => {
                    e.stopPropagation();
                    navigate(`/task/${task.task_id}`);
                  }}
                >
                  <span className="text-gray-300">
                    {task.title}
                    <span className="ml-2 text-xs text-gray-400">
                      {task.created_at &&
                        new Date(task.created_at).toLocaleString(undefined, {
                          month: 'short', day: 'numeric', year: 'numeric', hour: 'numeric', minute: '2-digit', hour12: true
                        })}
                      {task.archived &&
                        <>
                          {" (archived)"}
                        </>
                      }
                    </span>
                  </span>
                  <button
                    onClick={(e) => {
                      e.stopPropagation();
                      unarchiveTask(task.task_id);
                    }}
                    className="ml-1 px-1 rounded bg-green-700 text-white hover:bg-green-800 text-xs flex items-center"
                    style={{ marginTop: 0, marginBottom: 0, height: 'auto', minHeight: 0, paddingTop: 0, paddingBottom: 0, background: 'none', border: 'none', outline: 'none', boxShadow: 'none', cursor: 'pointer' }}
                    title="Unarchive task"
                  >
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor" className="w-4 h-4 mr-1">
                      <path strokeLinecap="round" strokeLinejoin="round" d="M12 19V6m0 0l-6 6m6-6l6 6" />
                    </svg>
                    Unarchive
                  </button>
                </li>
              ))}
            </ul>
          )
        )}
      </div>
      <div className="flex-1 flex flex-col justify-end">
        <div className="w-full flex justify-center pb-8">
          <form
            className="flex flex-col items-center w-full"
            style={{ maxWidth: '800px' }}
            onSubmit={e => {
              e.preventDefault();
              if (selectedEnv === "__add__" || selectedEnv === "__manage__") return;
              handleFormSubmit(e);
            }}
          >
            <div className="relative w-full" style={{ maxWidth: '800px', width: '100%' }}>
              <textarea
                rows={3}
                placeholder="Ask me your toxicology question [ Is green tea nephrotoxic? ]"
                value={chatInput}
                onChange={e => setChatInput(e.target.value)}
                onKeyDown={e => {
                  if (e.key === 'Enter' && !e.shiftKey) {
                    e.preventDefault();
                    if (selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
                      handleFormSubmit(e as any);
                    }
                  }
                }}
                className="w-full pt-4 pb-4 pr-28 pl-8 text-lg rounded-2xl border border-gray-700 bg-gray-900 bg-opacity-70 text-white resize-none min-h-[80px] shadow-2xl focus:outline-none focus:ring-2 focus:ring-green-400 text-left placeholder:text-left"
                style={{ minHeight: 80, fontFamily: 'inherit', width: '100%', boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)' }}
                disabled={selectedEnv === "__add__" || selectedEnv === "__manage__"}
              />
              <div className="absolute right-4 bottom-11 flex items-center space-x-2 z-30">
                <button
                  className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                  style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                  title="Upload CSV file"
                  onClick={handlePlusClick}
                  disabled={uploading}
                  type="button"
                >
                  <FaPlus className="w-5 h-5 text-black" />
                </button>
                <button
                  type="submit"
                  className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                  style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                  title="Submit"
                  disabled={selectedEnv === "__add__" || selectedEnv === "__manage__"}
                >
                  <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={2.2} stroke="black" className="w-5 h-5">
                    <path strokeLinecap="round" strokeLinejoin="round" d="M5 10l7-7m0 0l7 7m-7-7v18" />
                  </svg>
                </button>
              </div>
              <div className="absolute left-4 z-20 flex items-end" style={{ position: 'relative', width: '210px', overflow: 'visible', bottom: '3rem' }}>
                <div className="relative group" style={{ width: '100%', minWidth: '180px', maxWidth: '210px' }}>
                  <div className="flex items-center w-full" style={{ position: 'relative', height: '1.7rem', display: 'flex', alignItems: 'center' }}>
                    {/* Hidden span to measure width */}
                    <span
                      ref={spanRef}
                      style={{
                        position: "absolute",
                        visibility: "hidden",
                        whiteSpace: "nowrap",
                        fontWeight: "bold",
                        fontSize: "0.875rem",
                        fontFamily: "inherit",
                        padding: "0 16px"
                      }}
                    >
                      {selectedEnv === "__add__"
                        ? "+ Add environment"
                        : selectedEnv === "__manage__"
                          ? "⚙ Manage environments"
                          : (environments ?? []).find(e => e.environment_id === selectedEnv)
                            ? `env - ${(environments ?? []).find(e => e.environment_id === selectedEnv)?.title}`
                            : ""}
                    </span>
                    <select
                      ref={selectRef}
                      className={`font-bold text-white text-sm px-1 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-black group-hover:bg-opacity-60 group-hover:border group-hover:border-gray-700 group-hover:px-4 group-hover:pr-10 group-hover:cursor-pointer focus:bg-black focus:bg-opacity-60 focus:border focus:border-gray-700 focus:px-4 focus:pr-10 w-full ${loadingEnvironments ? 'opacity-50' : ''}`}
                      style={{
                        width: `${selectWidth}px`,
                        minWidth: "50px",
                        maxWidth: "260px",
                        height: '1.7rem',
                        borderRadius: '999px',
                        position: 'relative',
                        zIndex: 1,
                        backgroundColor: 'transparent',
                        cursor: 'pointer',
                        paddingRight: '1.5rem',
                      }}
                      value={selectedEnv}
                      onChange={handleDropdownChange}
                      disabled={loadingEnvironments}
                    >
                      {(environments ?? []).length > 0 && (environments ?? []).map(env => (
                        <option key={env.environment_id} value={env.environment_id} style={{ paddingLeft: '1rem' }}>
                          {`env - ${env.title}`}
                        </option>
                      ))}
                      <option value="__add__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>
                        {loadingEnvironments ? '+ Add environment' : '+ Add environment'}
                      </option>
                      <option value="__manage__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>&#9881; Manage environments</option>
                    </select>
                  </div>
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
      )}
    </div>
  );
};

export default Dashboard;
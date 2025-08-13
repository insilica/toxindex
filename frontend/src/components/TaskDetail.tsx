import React, { useEffect, useState, useCallback } from "react";
import { useParams, useNavigate } from "react-router-dom";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import HomeButton from "./shared/HomeButton";
import { FaDownload, FaFile, FaFilePdf, FaFileImage, FaFileExcel, FaFileWord, FaFileAlt, FaFileArchive, FaFileCode, FaFileCsv } from 'react-icons/fa';
import { MdOutlineTableChart } from 'react-icons/md';
import FilePreviewModal, { FilePreviewInline } from './shared/FilePreviewModal';
import { getWorkflowLabelById } from './shared/workflows';

interface Task {
  task_id: string;
  title: string;
  archived?: boolean;
  environment_id?: string;
  created_at?: string;
  finished_at?: string;
  session_id?: string;
  user_id?: string;
  workflow_id?: number;
}

const getDuration = (start?: string, end?: string) => {
  if (!start || !end) return null;
  const startTime = new Date(start).getTime();
  const endTime = new Date(end).getTime();
  const diff = Math.max(0, endTime - startTime);
  const mins = Math.floor(diff / 60000);
  const secs = Math.floor((diff % 60000) / 1000);
  return `${mins}m ${secs}s`;
};

// Helper to get icon by extension
const getFileIcon = (filename: string) => {
  const ext = filename.split('.').pop()?.toLowerCase();
  if (!ext) return <FaFile className="inline-block mr-1 w-4 h-4 align-text-bottom text-gray-400" />;

  if (["pdf"].includes(ext)) {
    return <FaFilePdf className="inline-block mr-1 w-4 h-4 align-text-bottom text-red-400" />;
  }
  if (["jpg", "jpeg", "png", "gif", "bmp", "svg"].includes(ext)) {
    return <FaFileImage className="inline-block mr-1 w-4 h-4 align-text-bottom text-orange-300" />;
  }
  if (["parquet"].includes(ext)) {
    return <MdOutlineTableChart className="inline-block mr-1 w-4 h-4 align-text-bottom text-gray-200" />;
  }
  if (["xls", "xlsx"].includes(ext)) {
    return <FaFileExcel className="inline-block mr-1 w-4 h-4 align-text-bottom text-green-400" />;
  }
  if (["csv"].includes(ext)) {
    return <FaFileCsv className="inline-block mr-1 w-4 h-4 align-text-bottom text-green-400" />;
  }
  if (["doc", "docx"].includes(ext)) {
    return <FaFileWord className="inline-block mr-1 w-4 h-4 align-text-bottom text-blue-400" />;
  }
  if (["txt"].includes(ext)) {
    return <FaFileAlt className="inline-block mr-1 w-4 h-4 align-text-bottom text-gray-300" />;
  }
  if (["zip", "rar", "tar", "gz"].includes(ext)) {
    return <FaFileArchive className="inline-block mr-1 w-4 h-4 align-text-bottom text-orange-400" />;
  }
  if (["js", "ts", "py", "json", "sh", "java", "cpp", "c", "cs", "rb"].includes(ext)) {
    return <FaFileCode className="inline-block mr-1 w-4 h-4 align-text-bottom text-purple-400" />;
  }
  return <FaFile className="inline-block mr-1 w-4 h-4 align-text-bottom text-gray-400" />;
};

const TaskDetail: React.FC = () => {
  const { task_id } = useParams<{ task_id: string }>();
  const [task, setTask] = useState<Task | null>(null);
  const [loading, setLoading] = useState(true);
  const [assistantMessage, setAssistantMessage] = useState<string | null>(null);
  const [messageLoading, setMessageLoading] = useState(false);
  const navigate = useNavigate();

  // Tab state
  // const [activeTab, setActiveTab] = useState<'MDrender' | 'MDraw' | 'JsonSchema'>('MDrender');
  // const [toxicitySchema, setToxicitySchema] = useState<any>(null);

  const [taskFiles, setTaskFiles] = useState<any[]>([]);
  const [envFiles, setEnvFiles] = useState<any[]>([]);

  const [previewFileId, setPreviewFileId] = useState<string | null>(null);
  const [previewOpen, setPreviewOpen] = useState(false);

  const fetchTask = useCallback(async () => {
    if (!task_id) return;
    setLoading(true);
    try {
      const res = await fetch(`/api/tasks/${task_id}`, { credentials: "include" });
      if (!res.ok) {
        setTask(null);
        setLoading(false);
        return;
      }
      const data = await res.json();
      if (data.error) {
        setTask(null);
      } else {
        setTask(data);
      }
    } catch (e) {
      setTask(null);
    } finally {
      setLoading(false);
    }
  }, [task_id]);

  // Fetch assistant message for this task
  useEffect(() => {
    const fetchAssistantMessage = async () => {
      if (!task_id) return;
      setMessageLoading(true);
      try {
        const res = await fetch(`/api/tasks/${task_id}/messages`, { credentials: "include" });
        if (!res.ok) throw new Error('Failed to fetch messages');
        const data = await res.json();
        if (data.messages && Array.isArray(data.messages)) {
          const assistantMessages = data.messages.filter((m: any) => m.role === "assistant");
          if (assistantMessages.length > 0) {
            setAssistantMessage(assistantMessages[assistantMessages.length - 1].content);
            return;
          }
        }
        setAssistantMessage(null); // No assistant message found
      } catch (e) {
        setAssistantMessage('Failed to load assistant message.');
        if (e instanceof Error) {
          console.error('Error loading assistant message:', e.message);
        } else {
          console.error('Unknown error loading assistant message:', e);
        }
      } finally {
        setMessageLoading(false);
      }
    };
    fetchAssistantMessage();
  }, [task_id]);

  // Fetch schema when JsonSchema tab is selected
  // useEffect(() => {
  //   if (activeTab === 'JsonSchema' && !toxicitySchema) {
  //     fetch('/api/schema/toxicity')
  //       .then(res => res.json())
  //       .then(setToxicitySchema)
  //       .catch(() => setToxicitySchema({ error: "Failed to load schema" }));
  //   }
  // }, [activeTab, toxicitySchema]);

  useEffect(() => {
    fetchTask();
  }, [fetchTask]);

  useEffect(() => {
    if (!task_id) return;
    // Fetch files for the task
    fetch(`/api/tasks/${task_id}/files`, { credentials: "include" })
      .then(res => res.json())
      .then(data => setTaskFiles(data.files || []))
      .catch(() => setTaskFiles([]));
  }, [task_id]);

  useEffect(() => {
    if (task && task.environment_id) {
      fetch(`/api/environments/${task.environment_id}/files`, { credentials: "include" })
        .then(res => res.json())
        .then(data => setEnvFiles(data.files || []))
        .catch(() => setEnvFiles([]));
    }
  }, [task]);

  // Find the latest two files by created_at (assuming taskFiles is sorted DESC by created_at)
  const latestFiles = taskFiles.slice(0, 2);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!task) return <div className="text-white p-8">Task not found.</div>;

  const duration = getDuration(task.created_at, task.finished_at);

  return (
    <div className="px-4 pt-4 pb-4 bg-neutral-800 text-neutral-200 flex flex-col min-h-screen w-full overflow-x-hidden" style={{ minHeight: '100vh' }}>
      <div className="flex items-center gap-4 px-4 py-4 border-b-2 border-gray-600 mb-6">
        <HomeButton />
        <span className="text-neutral-600 text-xl font-light mx-2">|</span>
        <span className="text-lg font-semibold text-white truncate">{task.title}</span>
        <span className="text-neutral-400 ml-4 text-sm whitespace-nowrap">
          {task.created_at ? new Date(task.created_at).toLocaleString() : "Unknown date"}
        </span>
      </div>
      <div className="flex flex-row w-full min-h-[60vh] overflow-hidden">
        {/* Sidebar: File List */}
        <div className="w-64 pr-4 border-r border-neutral-800 flex-shrink-0 overflow-y-auto" style={{ minWidth: 200, maxWidth: 250 }}>
          <div className="mb-6">
            <div className="mb-2 font-semibold text-neutral-200">Output files for this Task:</div>
            {taskFiles.length === 0 ? (
              <div className="text-gray-400 mb-2">No files for this task.</div>
            ) : (
              <ul className="mb-6">
                {(() => {
                  const maxFilenameLength = 32;
                  const truncate = (name: string) => {
                    if (name.length <= maxFilenameLength) return name;
                    const startLen = maxFilenameLength - 10;
                    const extMatch = name.match(/(\.[^./\\]+)$/);
                    const ext = extMatch ? extMatch[1] : '';
                    return `${name.slice(0, startLen)}...${ext}`;
                  };
                  return taskFiles.map(f => (
                    <li key={f.file_id} className="flex items-center space-x-3 mb-0.5">
                      <button
                        className="text-neutral-200 hover:text-blue-300 cursor-pointer bg-transparent border-none p-0 m-0"
                        title={f.filename}
                        style={{
                          background: 'none',
                          fontSize: '14px',
                          lineHeight: '1',
                          padding: '1px',
                          border: 'none',
                          height: '20px',
                          minHeight: 'unset',
                          fontWeight: 'normal',
                        }}
                        onClick={() => {
                          setPreviewFileId(f.file_id);
                          setPreviewOpen(true);
                        }}
                      >
                        {getFileIcon(f.filename)}
                        {truncate(f.filename)}
                      </button>
                      {/* <a
                        href={`/api/tasks/${f.task_id}/files/${f.file_id}/download`}
                        className="ml-2 px-1 py-1 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50"
                        title={`Download ${f.filename}`}
                        download
                      >
                        <FaDownload className="w-4 h-4" />
                      </a> */}
                    </li>
                  ));
                })()}
              </ul>
            )}
            <div className="mb-2 font-semibold text-neutral-200">Uploaded files:</div>
            {envFiles.length === 0 ? (
              <div className="text-neutral-400">No files for this environment.</div>
            ) : (
              <ul>
                {(() => {
                  const maxFilenameLength = 32;
                  const truncate = (name: string) => {
                    if (name.length <= maxFilenameLength) return name;
                    const startLen = maxFilenameLength - 10;
                    const extMatch = name.match(/(\.[^./\\]+)$/);
                    const ext = extMatch ? extMatch[1] : '';
                    return `${name.slice(0, startLen)}...${ext}`;
                  };
                  return envFiles.map(f => (
                    <li key={f.file_id} className="flex items-center space-x-2 mb-0.5">
                      <button
                        className="text-neutral-200 hover:text-blue-300 cursor-pointer bg-transparent border-none p-0 m-0"
                        title={f.filename}
                        style={{
                          background: 'none',
                          fontSize: '14px',
                          lineHeight: '1',
                          padding: '0',
                          border: 'none',
                          height: '20px',
                          minHeight: 'unset',
                          fontWeight: 'normal',
                        }}
                        onClick={() => {
                          setPreviewFileId(f.file_id);
                          setPreviewOpen(true);
                        }}
                      >
                        {getFileIcon(f.filename)}
                        {truncate(f.filename)}
                      </button>
                      {/* <a
                        href={`/api/environments/${f.environment_id}/files/${f.file_id}/download`}
                        className="ml-2 px-0 py-1 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50"
                        title={`Download ${f.filename}`}
                        download
                      >
                        <FaDownload className="w-4 h-4" />
                      </a> */}
                    </li>
                  ));
                })()}
              </ul>
            )}
          </div>
        </div>
        {/* Main Content */}
        <div className="flex-1 flex flex-col items-center justify-center w-full pl-4 overflow-hidden">
          {/* Tab content */}
          <div className="flex-1 flex flex-col items-center justify-center w-full max-w-full">
              <div className="w-full bg-neutral-800 rounded-lg p-4 shadow-inner overflow-hidden">
                {messageLoading ? (
                  <div className="text-gray-400">Loading result...</div>
                ) : assistantMessage ? (
                  <>
                    {latestFiles.length > 0 && latestFiles.map((file, idx) => (
                      <div key={file.file_id} className="mb-6 flex flex-col gap-2">
                        <div className="flex items-center gap-2 flex-wrap">
                          <span className="text-neutral-100 font-semibold whitespace-nowrap">Latest file{latestFiles.length > 1 ? ` #${idx + 1}` : ''}:</span>
                          <span className="text-neutral-300 font-mono text-sm truncate flex-1 min-w-0">{file.filename}</span>
                          <a
                            href={`/api/tasks/${file.task_id}/files/${file.file_id}/download`}
                            className="ml-2 px-2 py-1 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50 flex-shrink-0"
                            title={`Download ${file.filename}`}
                            download
                          >
                            <FaDownload className="w-3 h-3" />
                          </a>
                        </div>
                        <div className="border-b border-gray-700 my-2"></div>
                        {file.filename.toLowerCase().endsWith('.md') ? (
                          <div
                            className="prose prose-invert w-full cursor-pointer hover:bg-blue-950/30 transition overflow-x-auto"
                            style={{ lineHeight: 2, maxHeight: 650, overflowY: 'auto' }}
                            onClick={() => {
                              if (task?.session_id) {
                                navigate(`/chat/session/${task.session_id}?env=${task.environment_id}`);
                              }
                            }}
                            title="Go to chat session"
                          >
                            <ReactMarkdown
                              remarkPlugins={[remarkGfm]}
                              components={{
                                td: ({node, ...props}) => <td className="markdown-table-cell" {...props} />
                              }}
                            >
                              {assistantMessage}
                            </ReactMarkdown>
                          </div>
                        ) : (
                          <FilePreviewInline fileId={file.file_id} />
                        )}
                      </div>
                    ))}
                    {latestFiles.length === 0 && (
                      <div
                        className="prose prose-invert w-full cursor-pointer hover:bg-blue-950/30 transition overflow-x-auto"
                        style={{ lineHeight: 2, maxHeight: 650, overflowY: 'auto' }}
                        onClick={() => {
                          if (task?.session_id) {
                            navigate(`/chat/session/${task.session_id}?env=${task.environment_id}`);
                          }
                        }}
                        title="Go to chat session"
                      >
                        <ReactMarkdown
                          remarkPlugins={[remarkGfm]}
                          components={{
                            td: ({node, ...props}) => <td className="markdown-table-cell" {...props} />
                          }}
                        >
                          {assistantMessage}
                        </ReactMarkdown>
                      </div>
                    )}
                  </>
                ) : (
                  <div className="text-gray-400">No assistant message found for this task.</div>
                )}
              </div>

          </div>
        </div>
      </div>
      <footer className="mt-6 pt-1 pb-1 px-4 border-t border-gray-700 bg-gray-950 rounded-b-lg shadow-inner">
        <div className="w-full mx-auto">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-0 text-base overflow-x-auto min-w-0 w-full">
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Tool:</span>
              <span>{getWorkflowLabelById(task.workflow_id ?? 0) ?? <span className="text-gray-400">Unknown</span>}</span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Executed by user:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.user_id ? (
                  <button
                    className="text-purple-400 underline hover:text-purple-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/user/${task.user_id}`)}
                    title="Go to user profile"
                    aria-label="Go to user profile"
                  >
                    {task.user_id}
                  </button>
                ) : (
                  <span className="text-gray-400">Unknown</span>
                )}
              </span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Created:</span>
              <span>{task.created_at ? new Date(task.created_at).toLocaleString() : <span className="text-gray-400">Unknown</span>}</span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Finished:</span>
              <span>{task.finished_at ? new Date(task.finished_at).toLocaleString() : <span className="text-gray-400">-</span>}</span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Duration:</span>
              <span>{duration ?? <span className="text-gray-400">-</span>}</span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Environment:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.environment_id ? (
                  <button
                    className="text-blue-400 underline hover:text-blue-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/environments/details?env=${task.environment_id}`)}
                    title="Go to environment"
                    aria-label="Go to environment"
                  >
                    {task.environment_id}
                  </button>
                ) : (
                  <span className="text-gray-400">None</span>
                )}
              </span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Chat Session:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.session_id ? (
                  <button
                    className="text-green-400 underline hover:text-green-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/chat/session/${task.session_id}?env=${task.environment_id}`)}
                    title="Go to chat session"
                    aria-label="Go to chat session"
                  >
                    {task.session_id}
                  </button>
                ) : (
                  <span className="text-gray-400">None</span>
                )}
              </span>
            </div>
          </div>
        </div>
      </footer>
      <FilePreviewModal
        fileId={previewFileId}
        isOpen={previewOpen}
        onRequestClose={() => setPreviewOpen(false)}
      />
    </div>
  );
};

export default TaskDetail; 